#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int print_info = 0;

struct reduce_d_i {
    double max_d;
    int p_id;
} in_max, out_max;

struct reduce_i_i {
    int idx0;
    int p_id;
} in_min, out_min;

typedef struct tree_node {
    int center_idx;
    double radius;
    long left;
    long right;
    long node_id;
} node;

typedef struct pt_struct {
    int i;
    double proj;
    double *pt;
} pt;

int n_dims;
long np;
int p;
int p_initial;
int id;
int id_initial;
int new_id;
int aux;
double *p_aux;
int *block_size;
int *branch_size;
long *idx_global;
double *pt_arr_idx0;
double *pt_arr_a;
double *pt_arr_b;
double *proj;
double *pivots;
double *glob_pivots;
int *send_displs;
int *recv_displs;
int *send_counts;
int *recv_counts;
double *recv_buffer;
long *recv_buffer_long;
double *p_aux_2;
double *center;
double **centers;
long n_center = 0;
long tree_counter = 0;
int *sizes;
node *tree;
pt *pts;
int n_levels;
long id_last;

#define RANGE 10
#define BLOCK_LOW(id, p, np) ((id) * (np) / (p))
#define BLOCK_HIGH(id, p, np) (BLOCK_LOW((id) + 1, p, np) - 1)
#define BLOCK_SIZE(id, p, np) (BLOCK_HIGH(id, p, np) - BLOCK_LOW(id, p, np) + 1)

void print_point(double *pt, int dim) {
    for (int i = 0; i < dim; ++i) {
        fprintf(stdout, "%f ", pt[i]);
        fflush(stdout);
    }
    printf("\n");
    fflush(stdout);
}

void swap(int i, int j) {
    double temp = proj[i];
    proj[i] = proj[j];
    proj[j] = temp;
    double temp_pt;
    for (int k = 0; k < n_dims; k++) {
        temp_pt = p_aux[i * n_dims + k];
        p_aux[i * n_dims + k] = p_aux[j * n_dims + k];
        p_aux[j * n_dims + k] = temp_pt;
    }
    long idx_temp = idx_global[i];
    idx_global[i] = idx_global[j];
    idx_global[j] = idx_temp;
}

double dist(double *pt1, double *pt2) {
    double dist = 0.0;
    for (int d = 0; d < n_dims; ++d)
        dist += (pt1[d] - pt2[d]) * (pt1[d] - pt2[d]);
    return dist;
}

void swap_vec(int i, int j) {
    double temp = glob_pivots[i];
    glob_pivots[i] = glob_pivots[j];
    glob_pivots[j] = temp;
}

int partition(double pivot, int l, int r) {

    if (r == l + 1) {
        if (proj[l] <= pivot)
            return l + 1;
        else
            return l;
    }

    int i = l - 1;
    int j;
    for (j = l; j < r; j++) {
        if (proj[j] <= pivot) {
            i++;
            swap(i, j);
        }
    }

    return i + 1;
}

int is_seq(int i, int j) {
    return proj[i] <= proj[j];
}

int is_seq_2(int i, int j) {
    return glob_pivots[i] <= glob_pivots[j];
}

double quickselect_seq_2(int k, int l, int h) {

    if (h == l + 1)
        return glob_pivots[l];

    else {
        int i = l;
        int j = h;

        while (i < j) {
            do {
                i++;
            } while (i < h && is_seq_2(i, l));
            if (i == h) {
                j = h - 1;
                break;
            }
            do {
                j--;
            } while (!is_seq_2(j, l));
            if (i < j)
                swap_vec(i, j);
        }

        swap_vec(l, j);

        if (j - l == k)
            return glob_pivots[j];
        if (k < j - l)
            return quickselect_seq_2(k, l, j);
        else
            return quickselect_seq_2(k - (j - l) - 1, j + 1, h);
    }
}

int quickselect_seq(int k, int l, int h) {

    if (h == l + 1)
        return l;

    else {
        int i = l;
        int j = h;

        while (i < j) {
            do {
                i++;
            } while (i < h && is_seq(i, l));
            if (i == h) {
                j = h - 1;
                break;
            }
            do {
                j--;
            } while (!is_seq(j, l));
            if (i < j)
                swap(i, j);
        }

        swap(l, j);

        if (j - l == k)
            return j;
        if (k < j - l)
            return quickselect_seq(k, l, j);
        else
            return quickselect_seq(k - (j - l) - 1, j + 1, h);
    }
}

double orth_projv1(double *a, double *b, double *pt) {
    double proj = 0.0;
    for (int i = 0; i < n_dims; ++i)
        proj += (pt[i] - a[i]) * (b[i] - a[i]);
    return proj;
}

void swap_seq(long *i, long *j) {
    long temp = *i;
    *i = *j;
    *j = temp;
}

void furthest_apart(int l, int r, long *a_idx, long *b_idx) {
    double max_d = 0.0;
    double d;
    long idx0 = l;
    long idx0_comp = pts[l].i;
    for (int i = l + 1; i < r; ++i) {
        if (pts[i].i < idx0_comp) {
            idx0 = i;
            idx0_comp = pts[i].i;
        }
    }
    for (int i = l; i < r; ++i) {
        d = dist(pts[idx0].pt, pts[i].pt);
        if (d > max_d) {
            max_d = d;
            *a_idx = i;
        }
    }
    max_d = 0;
    for (int i = l; i < r; ++i) {
        d = dist(pts[*a_idx].pt, pts[i].pt);
        if (d > max_d) {
            max_d = d;
            *b_idx = i;
        }
    }
    if (pts[*a_idx].pt[0] > pts[*b_idx].pt[0])
        swap_seq(a_idx, b_idx);
}

int is_seq_seq(int i, int j) {
    return pts[i].proj <= pts[j].proj;
}

void swap_pt(long i, long j) {
    pt temp = pts[i];
    pts[i] = pts[j];
    pts[j] = temp;
}

int get_kth_element(int k, int l, int h) {

    if (h == l + 1)
        return l;

    else {
        int i = l;
        int j = h;

        while (i < j) {
            do {
                i++;
            } while (i < h && is_seq_seq(i, l));
            if (i == h) {
                j = h - 1;
                break;
            }
            do {
                j--;
            } while (!is_seq_seq(j, l));
            if (i < j)
                swap_pt(i, j);
        }

        swap_pt(l, j);

        if (j - l == k)
            return j;
        if (k < j - l)
            return get_kth_element(k, l, j);
        else
            return get_kth_element(k - (j - l) - 1, j + 1, h);
    }
}

void get_median(int l, int h, int *median_1, int *median_2) {
    if ((h - l) % 2 == 1)
        *median_1 = get_kth_element((h - l - 1) / 2, l, h);
    else {
        *median_1 = get_kth_element((h - l) / 2 - 1, l, h);
        *median_2 = l + (h - l) / 2;
        for (int i = l + (h - l) / 2 + 1; i < h; ++i)
            if (is_seq_seq(i, *median_2))
                *median_2 = i;
    }
}

// long tree_idx;

// This is the algorithm sequential version but with tasks in the recursive calls
void ballAlg_tasks(long l, long r, long n_id, int lvl) {

    long temp_counter;

    if (r - l == 1) {
        #pragma omp critical
        {
            temp_counter = tree_counter++;
        }
        tree[temp_counter].center_idx = l;
        tree[temp_counter].node_id = n_id;
        tree[temp_counter].radius = 0;
        tree[temp_counter].left = -1;
        tree[temp_counter].right = -1;
        return;
    }

    long a, b;
    long c_id = -1;
    // 2. Compute points a and b, furthest apart in the current set (approx)
    furthest_apart(l, r, &a, &b);

    // 3. Perform the orthogonal projection of all points onto line ab
    for (int i = l; i < r; ++i)
        pts[i].proj = orth_projv1(pts[a].pt, pts[b].pt, pts[i].pt);

    // 4. Compute the center, defined as the median point over all projections
    int m1, m2 = -1;
    double u, aux;
    double abnorm = 0;
    get_median(l, r, &m1, &m2);
    #pragma omp critical
    {
        temp_counter = tree_counter++;
        c_id = n_center++;
    }

    if ((r - l) % 2)
        u = pts[m1].proj;
    else
        u = (pts[m1].proj + pts[m2].proj) / 2;

    for (int i = 0; i < n_dims; ++i) {
        aux = (pts[b].pt[i] - pts[a].pt[i]);
        centers[c_id][i] = u * aux;
        abnorm += aux * aux;
    }
    for (int i = 0; i < n_dims; ++i)
        centers[c_id][i] = centers[c_id][i] / abnorm + pts[a].pt[i];

    double max_r = 0, rad;
    for (int i = l; i < r; ++i) {
        rad = dist(centers[c_id], pts[i].pt);
        if (rad > max_r)
            max_r = rad;
    }
    tree[temp_counter].radius = sqrt(max_r);
    tree[temp_counter].center_idx = c_id;
    tree[temp_counter].node_id = n_id;

    if (((block_size[id_initial] & (block_size[id_initial] - 1)) != 0) && lvl == n_levels - 1) {
        #pragma omp critical
        {
            tree[temp_counter].left = id_last;
            tree[temp_counter].right = id_last + 1;
            id_last += 2;
        }
    } else {
        tree[temp_counter].left = 2 * n_id + 1;
        tree[temp_counter].right = 2 * n_id + 2;
    }

    #pragma omp task
        ballAlg_tasks(l, l + (r - l) / 2, tree[temp_counter].left, lvl + 1);
    #pragma omp task
        ballAlg_tasks(l + (r - l) / 2, r, tree[temp_counter].right, lvl + 1);
}

// Sequential version of the algorithm with no tasks in the recursive calls
void ballAlg(long l, long r, long n_id, int lvl) {

    long temp_counter;

    if (r - l == 1) {
        #pragma omp critical
        {
            temp_counter = tree_counter++;
        }
        tree[temp_counter].center_idx = l;
        tree[temp_counter].radius = 0;
        tree[temp_counter].node_id = n_id;
        tree[temp_counter].left = -1;
        tree[temp_counter].right = -1;
        return;
    }

    long a, b;
    long c_id = -1;
    // 2. Compute points a and b, furthest apart in the current set (approx)
    furthest_apart(l, r, &a, &b);
    double * pt_a = pts[a].pt;
    double * pt_b = pts[b].pt;

    // 3. Perform the orthogonal projection of all points onto line ab
    for (int i = l; i < r; ++i)
        pts[i].proj = orth_projv1(pt_a, pt_b, pts[i].pt);

    // 4. Compute the center, defined as the median point over all projections
    int m1, m2 = -1;
    double u, aux;
    double abnorm = 0;
    get_median(l, r, &m1, &m2);
    #pragma omp critical
    {
        c_id = n_center++;
        temp_counter = tree_counter++;
    }

    if ((r - l) % 2)
        u = pts[m1].proj;
    else
        u = (pts[m1].proj + pts[m2].proj) / 2;

    for (int i = 0; i < n_dims; ++i) {
        aux = (pt_b[i] - pt_a[i]);
        centers[c_id][i] = u * aux;
        abnorm += aux * aux;
    }
    for (int i = 0; i < n_dims; ++i)
        centers[c_id][i] = centers[c_id][i] / abnorm + pt_a[i];

    double max_r = 0, rad;
    for (int i = l; i < r; ++i) {
        rad = dist(centers[c_id], pts[i].pt);
        if (rad > max_r)
            max_r = rad;
    }
    tree[temp_counter].radius = sqrt(max_r);
    tree[temp_counter].center_idx = c_id;
    tree[temp_counter].node_id = n_id;

    if (((block_size[id_initial] & (block_size[id_initial] - 1)) != 0) && lvl == n_levels - 1) {
        // printf("id: %d, node_id: %ld, id_last: %ld\n", id_initial, n_id, id_last);
        // fflush(stdout);
        #pragma omp critical
        {
            tree[temp_counter].left = id_last;
            tree[temp_counter].right = id_last + 1;
            id_last += 2;
        }
    } else {
        tree[temp_counter].left = 2 * n_id + 1;
        tree[temp_counter].right = 2 * n_id + 2;
    }

    ballAlg(l, l + (r - l) / 2, tree[temp_counter].left, lvl + 1);
    ballAlg(l + (r - l) / 2, r, tree[temp_counter].right, lvl + 1);
}

int max_parallel_level; // Max level until which paralelization inside the algorithm will be used

/* Parallel algoritm to be used in the first levels */
void ballAlg_omp(long l, long r, long n_id, long lvl, int threads) {

    // printf("id: %d, id_last: %ld\n", id_initial, id_last);
    // fflush(stdout);

    long a = -1;         // Global a
    long b = -1;         // Global b
    int m1, m2 = -1;     // Global median indexes
    double abnorm = 0.0; // Global variable to be used in the reduction for

    long idx0, idx0_comp = np;
    double max_d = 0.0; // Global max distance
    double max_r = 0.0; // Global max radius
    long c_id = -1;     // Index for the centers vector
    double u;
    long temp_counter; // Index for the current tree struct position

    #pragma omp parallel num_threads(threads) private(pt_arr_a, pt_arr_b)
    {

        if (r - l == 1) {
            #pragma omp single
            {
                #pragma omp critical
                temp_counter = tree_counter++;

                tree[temp_counter].center_idx = l;
                tree[temp_counter].node_id = n_id;
                tree[temp_counter].radius = 0;
                tree[temp_counter].left = -1;
                tree[temp_counter].right = -1;
            }
        }

        else {
            long idx0_t = l;             // Minimum id inside each thread
            long idx0_t_glob = pts[l].i; // Minimum id inside each thread
            #pragma omp for
            for (int i = l + 1; i < r; ++i) {
                if (pts[i].i < idx0_t_glob) {
                    idx0_t = i;
                    idx0_t_glob = pts[i].i;
                }
            }

            #pragma omp critical // Calculation of the minimum of minimums
            {
                if (idx0_t_glob < idx0_comp) {
                    idx0 = idx0_t; // Minimum id from all threads
                    idx0_comp = idx0_t_glob;
                }
            }

            #pragma omp barrier // Do not let any threads move foward before we complete the calculation of idx0, since it's going to be used next

            double d;
            double max_d_t = 0.0; // Maximum distance for each thread
            long a_t = -1;        // Corresponding index for each thread
            #pragma omp for
            for (int i = l; i < r; ++i) {
                d = dist(pts[idx0].pt, pts[i].pt);
                if (d > max_d_t) {
                    max_d_t = d;
                    a_t = i;
                }
            }

            #pragma omp critical // The global maximum is the maximum of the threads maximums
            {
                if (max_d_t > max_d) {
                    max_d = max_d_t;
                    a = a_t; // Also assign the global id
                }
            }

            #pragma omp barrier
            max_d = 0; // Reset the global distance

            long b_t = -1;
            max_d_t = 0; // Reset the private distances
            #pragma omp for
            for (int i = l; i < r; ++i) {
                d = dist(pts[a].pt, pts[i].pt);
                if (d > max_d_t) {
                    max_d_t = d;
                    b_t = i;
                }
            }

            #pragma omp critical
            {
                if (max_d_t > max_d) {
                    max_d = max_d_t;
                    b = b_t;
                }
            }
            #pragma omp barrier

            #pragma omp single
            {
                if (pts[a].pt[0] > pts[b].pt[0])
                    swap_seq(&a, &b);
                // printf("id: %ld, pts[a].pt: %f %f %f, pts[b].pt: %f %f %f\n", n_id, pts[a].pt[0], pts[a].pt[1], pts[a].pt[2], pts[b].pt[0], pts[b].pt[1], pts[b].pt[2]);
                // fflush(stdout);
            }
            pt_arr_a = pts[a].pt;
            pt_arr_b = pts[b].pt;

// 3. Perform the orthogonal projection of all points onto line ab
            #pragma omp for
            for (int i = l; i < r; ++i) {
                pts[i].proj = orth_projv1(pts[a].pt, pts[b].pt, pts[i].pt);
            }

// 4. Compute the center, defined as the median point over all projections
            #pragma omp single
            {
                get_median(l, r, &m1, &m2);
                // if ((r - l) % 2)
                //     printf("id: %ld, pt[m1]: %f %f %f\n", n_id, pts[m1].pt[0], pts[m1].pt[1], pts[m1].pt[2]);
                // else
                //     printf("id: %ld, pt[m1]: %f %f %f, pt[m2]: %f %f %f\n", n_id, pts[m1].pt[0], pts[m1].pt[1], pts[m1].pt[2], pts[m2].pt[0], pts[m2].pt[1], pts[m2].pt[2]);
                #pragma omp critical // Increment the counter of centers
                {
                    c_id = n_center++;
                    temp_counter = tree_counter++;
                }
            }

            #pragma omp single
            {
                if ((r - l) % 2)
                    u = pts[m1].proj;
                else {
                    u = (pts[m1].proj + pts[m2].proj) / 2;
                }
                // printf("id: %ld, u: %f\n", n_id, u);
                // fflush(stdout);
            }

            // printf("id: %ld, pt_arr_a: %f %f %f, pt_arr_b: %f %f %f\n", n_id, pt_arr_a[0], pt_arr_a[1], pt_arr_a[2], pt_arr_b[0], pt_arr_b[1], pt_arr_b[2]);
            // fflush(stdout);
            double aux;
            #pragma omp for reduction(+ : abnorm) // Using reduction since abnorm is a global variable
            for (int i = 0; i < n_dims; ++i) {
                aux = (pt_arr_b[i] - pt_arr_a[i]);
                centers[c_id][i] = u * aux;
                abnorm += aux * aux;
            }

            #pragma omp for
            for (int i = 0; i < n_dims; ++i)
                centers[c_id][i] = centers[c_id][i] / abnorm + pt_arr_a[i];

            // #pragma omp single
            // {
            //     printf("id: %ld, abnorm: %f, centers[c_id]: %f %f %f\n", n_id, abnorm, centers[c_id][0], centers[c_id][1], centers[c_id][2]);
            //     fflush(stdout);
            // }

            double rad;
            double max_r_t = 0; // Maximum radius inside each thread

            #pragma omp for
            for (int i = l; i < r; ++i) {
                rad = dist(centers[c_id], pts[i].pt);
                if (rad > max_r_t)
                    max_r_t = rad;
            }

            #pragma omp critical // Get the maximum radius from the radius in each thread
            {
                if (max_r_t > max_r)
                    max_r = max_r_t;
            }
            #pragma omp barrier

            #pragma omp single
            {
                tree[temp_counter].radius = sqrt(max_r);
                tree[temp_counter].center_idx = c_id;
                tree[temp_counter].node_id = n_id;
                // If the number of points is not a power of two and if we are in the last level of the tree, use variavle id_last to calculate the id's
                if (((block_size[id_initial] & (block_size[id_initial] - 1)) != 0) && lvl == n_levels - 1) {
                    #pragma omp critical
                    {
                        tree[temp_counter].left = id_last;
                        tree[temp_counter].right = id_last + 1;
                        id_last += 2;
                    }
                }

                else {
                    tree[temp_counter].left = 2 * n_id + 1;
                    tree[temp_counter].right = 2 * n_id + 2;
                }

                // If we are in the last level with paralelization
                if (lvl == max_parallel_level) {
                    // Call sequential version with tasks if the number of threads is not a power of 2
                    if (lvl == 0 && (threads & (threads - 1)) != 0) {
                        #pragma omp task
                            ballAlg_tasks(l, l + (r - l) / 2, tree[temp_counter].left, lvl + 1);
                        #pragma omp task
                            ballAlg_tasks(l + (r - l) / 2, r, tree[temp_counter].right, lvl + 1);
                    }
                    // Call sequential version with no tasks
                    else {
                        #pragma omp task
                            ballAlg(l, l + (r - l) / 2, tree[temp_counter].left, lvl + 1);
                        #pragma omp task
                            ballAlg(l + (r - l) / 2, r, tree[temp_counter].right, lvl + 1);
                    }
                } else {
                    #pragma omp task
                        ballAlg_omp(l, l + (r - l) / 2, tree[temp_counter].left, lvl + 1, threads / 2);
                    #pragma omp task
                        ballAlg_omp(l + (r - l) / 2, r, tree[temp_counter].right, lvl + 1, threads - threads / 2);
                }
            }
        }
    }
}

/*
long ballAlg_seq(long l, long r) {

    long id_temp = tree_idx++;

    if (r - l == 1) {
        tree[tree_counter].center_idx = l;
        tree[tree_counter].radius = 0;
        tree[tree_counter].left = -1;
        tree[tree_counter].right = -1;
        tree[tree_counter].node_id = id_temp;
        tree_counter++;
        return id_temp;
    }

    long temp_counter = tree_counter++;

    long a, b;
    furthest_apart(l, r, &a, &b);

    pt_arr_a = pts[a].pt;
    pt_arr_b = pts[b].pt;

    for (int i = l; i < r; ++i)
        pts[i].proj = orth_projv1(pt_arr_a, pt_arr_b, pts[i].pt);

    int m1, m2 = -1;
    get_median(l, r, &m1, &m2);

    double abnorm = 0.0, aux, u;
    if ((r - l) % 2)
        u = pts[m1].proj;
    else
        u = (pts[m1].proj + pts[m2].proj) / 2;

    for (int i = 0; i < n_dims; ++i) {
        aux = (pt_arr_b[i] - pt_arr_a[i]);
        centers[n_center][i] = u * aux;
        abnorm += aux * aux;
    }
    for (int i = 0; i < n_dims; ++i)
        centers[n_center][i] = centers[n_center][i] / abnorm + pt_arr_a[i];

    double max_r = 0, rad;
    for (int i = l; i < r; ++i) {
        rad = dist(centers[n_center], pts[i].pt);
        if (rad > max_r)
            max_r = rad;
    }

    tree[temp_counter].center_idx = n_center++;
    tree[temp_counter].radius = sqrt(max_r);
    tree[temp_counter].node_id = id_temp;

    tree[temp_counter].left = ballAlg_seq(l, l + (r - l) / 2);
    tree[temp_counter].right = ballAlg_seq(l + (r - l) / 2, r);

    return id_temp;
}
*/

void ballAlg_mpi(long n_points, long n_id, int lvl, MPI_Comm comm, int threads) {

    // Only called if no leaves will be done here

    if (print_info == 1) {
        printf("block_size[id_initial]: %d, np: %ld, id_inital: %d, id: %d, lvl: %d, n_id: %ld\n", block_size[id_initial], n_points, id_initial, id, lvl, n_id);
        fflush(stdout);
    }

    /* LOAD BALACING: EVENS OUT THE NUMBER OF POINTS IN THE PROCESSORS WITH THE LARGEST AND SMALLEST NUMBER OF POINTS*/
    if (n_id != 0) {
        // printf("Before - block_size[id_initial]: %d, np: %ld, id_inital: %d, id: %d, lvl: %d, n_id: %ld\n", block_size[id_initial], n_points, id_initial, id, lvl, n_id);
        // fflush(stdout);

        MPI_Allgather(&block_size[id_initial], 1, MPI_INT, sizes, 1, MPI_INT, comm);

        int p_max = 0;
        int p_min = 0;
        int max = sizes[0];
        int min = sizes[0];
        for (int i = 1; i < p; i++) {
            if (sizes[i] > max) {
                max = sizes[i];
                p_max = i;
            }
            if (sizes[i] < min) {
                min = sizes[i];
                p_min = i;
            }
        }

        if (id == p_max) {
            int send_quant = (max - min) / 2;
            int first_id = block_size[id_initial] - send_quant;
            MPI_Send(proj + first_id, send_quant, MPI_DOUBLE, p_min, 0, comm);
            MPI_Send(p_aux + first_id * n_dims, (send_quant)*n_dims, MPI_DOUBLE, p_min, 1, comm);
            MPI_Send(idx_global + first_id, send_quant, MPI_LONG, p_min, 2, comm);
            block_size[id_initial] -= send_quant;
        }
        if (id == p_min) {
            int recv_quant = (max - min) / 2;
            int first_id = block_size[id_initial];
            MPI_Recv(proj + first_id, recv_quant, MPI_DOUBLE, p_max, 0, comm, MPI_STATUS_IGNORE);
            MPI_Recv(p_aux + first_id * n_dims, (recv_quant)*n_dims, MPI_DOUBLE, p_max, 1, comm, MPI_STATUS_IGNORE);
            MPI_Recv(idx_global + first_id, recv_quant, MPI_LONG, p_max, 2, comm, MPI_STATUS_IGNORE);
            block_size[id_initial] += recv_quant;
        }

        // printf("After - block_size[id_initial]: %d, np: %ld, id_inital: %d, id: %d, lvl: %d, n_id: %ld\n", block_size[id_initial], n_points, id_initial, id, lvl, n_id);
        // fflush(stdout);
    }

    /* FINDS THE POINT WITH THE LOWEST INITIAL INDEX (LOCALLY) */
    in_min.idx0 = idx_global[0];
    in_min.p_id = id;
    long local = 0;
    for (int i = 1; i < block_size[id_initial]; i++) {
        if (idx_global[i] < in_min.idx0) {
            in_min.idx0 = idx_global[i];
            local = i;
        }
    }
    // printf("id: %d, in_min.idx0: %d, in_min.p_id: %d, local: %d\n", id, in_min.idx0, in_min.p_id, local);
    // fflush(stdout);

    /* REDUCES OVER ALL THE PROCESSORS, GETTING THE GLOBAL MINIMUM */
    MPI_Allreduce(&in_min, &out_min, 1, MPI_2INT, MPI_MINLOC, comm);

    /* THE PROCESSOR WITH THE POINT WITH THE LOWEST INITIAL INDEX WILL SEND IT TO ALL THE OTHER PROCESSORS */
    if (id == out_min.p_id) {
        for (int j = 0; j < n_dims; j++) {
            pt_arr_idx0[j] = p_aux[local * n_dims + j];
        }
    }
    MPI_Bcast(pt_arr_idx0, n_dims, MPI_DOUBLE, out_min.p_id, comm);

    /* FINDS THE LOCAL POINT FURTHEST APART FROM THE POINT WITH LOWER INDEX */
    double d;
    int a = 0;
    in_max.max_d = 0;
    in_max.p_id = id;
    for (int i = 0; i < block_size[id_initial]; ++i) {
        d = 0;
        /* Inline distance calculation */
        for (int j = 0; j < n_dims; j++)
            d += (pt_arr_idx0[j] - p_aux[i * n_dims + j]) * (pt_arr_idx0[j] - p_aux[i * n_dims + j]);
        if (d > in_max.max_d) {
            in_max.max_d = d;
            a = i;
        }
    }
    /* FINDS THE GLOBAL POINT FURTHEST APART FROM THE POINT WITH LOWER INDEX - point a */
    MPI_Allreduce(&in_max, &out_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);

    /* THE PROCESSOR WITH THE a POINT SENDS IT TO ALL THE OTHER PROCESSORS */
    if (id == out_max.p_id) {
        for (int j = 0; j < n_dims; j++) {
            pt_arr_a[j] = p_aux[a * n_dims + j];
        }
    }
    MPI_Bcast(pt_arr_a, n_dims, MPI_DOUBLE, out_max.p_id, comm);

    /* SAME AS BEFORE, BUT FOR THE POINT b */
    int b = 0;
    in_max.max_d = 0;
    in_max.p_id = id;
    for (int i = 0; i < block_size[id_initial]; ++i) {
        d = 0;
        for (int j = 0; j < n_dims; j++)
            d += (pt_arr_a[j] - p_aux[i * n_dims + j]) * (pt_arr_a[j] - p_aux[i * n_dims + j]);
        if (d > in_max.max_d) {
            in_max.max_d = d;
            b = i;
        }
    }

    MPI_Allreduce(&in_max, &out_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);

    if (id == out_max.p_id) {
        for (int j = 0; j < n_dims; j++) {
            pt_arr_b[j] = p_aux[b * n_dims + j];
        }
    }
    MPI_Bcast(pt_arr_b, n_dims, MPI_DOUBLE, out_max.p_id, comm);

    /* ALL PROCESSORS NOW HAVE THE POINTS a AND b */

    /* SWAPS a AND b, IF NEEDED, SO THAT a HAS A LOWER FIRST COORDINATE THAN b */
    if (pt_arr_a[0] > pt_arr_b[0]) {
        double temp;
        for (int i = 0; i < n_dims; i++) {
            temp = pt_arr_a[i];
            pt_arr_a[i] = pt_arr_b[i];
            pt_arr_b[i] = temp;
        }
    }

    /* EACH PROCESSOR COMPUTES ITS PROJECTIONS OVER THE ab LINE */
    for (int i = 0; i < block_size[id_initial]; ++i) {
        proj[i] = 0;
        for (int j = 0; j < n_dims; j++)
            proj[i] += (p_aux[i * n_dims + j] - pt_arr_a[j]) * (pt_arr_b[j] - pt_arr_a[j]);
    }

    /* START OF THE PARALLEL MEDIAN FINDING, USING A VARIATION OF PARALLEL SORTING WITH REGULAR SAMPLING (PSRS) */

    /* FIRST PIVOT: POINT WITH THE SMALLER PROJECTION IN EACH PROCESSOR */
    int min_index = 0;
    for (int i = 0; i < block_size[id_initial]; ++i)
        if (proj[i] < proj[min_index])
            min_index = i;
    swap(0, min_index);
    pivots[0] = proj[0];

    int temp_idx = -1;
    int piv_idx, prev_piv_idx = 0;
    /* FINDS THE REMAINING p-1 PIVOTS (p = # of processors), SO THAT THEY ARE AS EVENLY SPACED AS POSSIBLE */
    for (int i = 1; i < p; ++i) {
        piv_idx = i * block_size[id_initial] / p - prev_piv_idx - 1;
        temp_idx = quickselect_seq(piv_idx, prev_piv_idx + 1, block_size[id_initial]);
        prev_piv_idx = temp_idx;
        pivots[i] = proj[temp_idx];
    }

    /* GATHERS ALL PIVOTS IN THE PROCESS WITH INITIAL RANK 0 */
    MPI_Gather(pivots, p, MPI_DOUBLE, glob_pivots, p, MPI_DOUBLE, 0, comm);

    /* PROCESS WITH INITIAL RANK 0 FINDS p-1 PIVOTS FROM THE p^2 PIVOTS, ALSO EVENLY SPACED */
    if (!id) {
        prev_piv_idx = -1;
        for (int i = 1; i < p; ++i) {
            piv_idx = i * p - prev_piv_idx - 1;
            pivots[i - 1] = quickselect_seq_2(piv_idx, prev_piv_idx + 1, p * p);
            prev_piv_idx += piv_idx + 1;
        }
    }

    /* PROCESS WITH INITIAL RANK 0 SENDS THE FOUND PIVOTS TO THE REMAINING PROCESSES */
    MPI_Bcast(pivots, p - 1, MPI_DOUBLE, 0, comm);

    /* EACH PROCESS PARTITIONS ITS DATA ACCORDING TO THE RECENTLY FOUND p-1 PIVOTS */
    /* THEY ALSO KEEP A COUNTER OF THE SIZE OF EACH BLOCK */
    send_counts[0] = partition(pivots[0], 0, block_size[id_initial]);
    int sum = send_counts[0];
    for (int i = 1; i < p - 1; i++) {
        send_counts[i] = partition(pivots[i], sum, block_size[id_initial]) - sum;
        sum += send_counts[i];
    }
    send_counts[p - 1] = block_size[id_initial] - sum;

    /* ALL DATA LOWER THAN THE FIRST PIVOT WILL BE SENT TO THE PROCESS WITH THE LOWEST RANK */
    /* DATA BETWEEN THE FIRST AND SECOND PIVOTS WILL BE SENT TO THE PROCESS WITH THE SECOND LOWEST RANK ... */
    /* ... DATA GREATER THAN THE LAST PIVOT WILL BE SENT TO THE PROCESS WITH HIGHEST RANK */

    /* EACH PROCESS STARTS BY SENDING THE SIZE OF EACH OF ITS CHUNKS */
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, comm);

    int sum1 = 0; //
    int sum2 = 0;
    for (int i = 0; i < p; i++) {
        send_displs[i] = sum2;  // where the data from this processor is partitioned to be sent to the remaining processes
        recv_displs[i] = sum1;  // where the data from the other processes will be placed in the current process
        sum1 += recv_counts[i]; // total receaving size
        sum2 += send_counts[i]; // total sending size
    }

    if (print_info == 1) {
        printf("id: %d, sum1: %d\n", id, sum1);
        fflush(stdout);
    }

    /* ALL PROCESSES EXCHANGE CHUNKS SO THAT ALL THE DATA IN PROCESS i IS GREATER THAN THAT OF i-1 AND SMALLER THAN THE DATA IN i+1 */
    /* PROJECTIONS EXCHANGE */
    MPI_Alltoallv(proj, send_counts, send_displs, MPI_DOUBLE, recv_buffer, recv_counts, recv_displs, MPI_DOUBLE, comm);
    double *temp = proj;
    proj = recv_buffer;
    recv_buffer = temp;

    /* GLOBAL INDECES EXCHANGE */
    MPI_Alltoallv(idx_global, send_counts, send_displs, MPI_LONG, recv_buffer_long, recv_counts, recv_displs, MPI_LONG, comm);
    long *temp_long = idx_global;
    idx_global = recv_buffer_long;
    recv_buffer_long = temp_long;

    /* ACTUAL POINTS EXCHANGE */
    for (int i = 0; i < p; ++i) {
        send_counts[i] *= n_dims;
        send_displs[i] *= n_dims;
        recv_counts[i] *= n_dims;
        recv_displs[i] *= n_dims;
    }
    MPI_Alltoallv(p_aux, send_counts, send_displs, MPI_DOUBLE, p_aux_2, recv_counts, recv_displs, MPI_DOUBLE, comm);
    temp = p_aux;
    p_aux = p_aux_2;
    p_aux_2 = temp;

    /* THE MEDIAN IS, WITH AN HIGH PROBABILITY, IN THE MIDDLE PROCESSOR (OR ONE OF THE MIDDLE, FOR EVEN NUMBER OF PROCESSORS) */

    /* SUMS THE SIZES FROM ALL THE LEFT PROCESSORS, UP TO p/2 (IF EVEN NUMBER OF PROCESSORS) OR p/2 + 1 (IF ODD) */
    int size = (id < p / 2 + (p % 2)) ? sum1 : 0;
    int total_size = 0;
    MPI_Allreduce(&size, &total_size, 1, MPI_INT, MPI_SUM, comm);

    if (print_info == 1) {
        printf("id: %d, size: %d, total_size: %d\n", id, size, total_size);
        fflush(stdout);
    }

    /* FINDS THE ACTUAL MEDIAN (FINALLY) */
    int m1 = -1, m2 = -1;
    double min_proj;
    double u;
    if ((id == p / 2) && (p % 2 == 1)) {
        if ((n_points / 2 - !(n_points % 2) < (total_size - sum1)) || (n_points / 2 - !(n_points % 2) >= total_size)) {
            fprintf(stderr, "MEDIAN OUT OF BOUNDS - m1: %ld, sum1: %d\n", n_points / 2 - !(n_points % 2) - (total_size - sum1), sum1);
            fflush(stderr);
        }
        m1 = quickselect_seq(n_points / 2 - !(n_points % 2) - (total_size - sum1), 0, sum1);
        if (n_points % 2 == 0) {
            min_proj = proj[m1 + 1];
            m2 = m1 + 1;
            for (int i = m1 + 1; i < sum1; i++)
                if (proj[i] < min_proj) {
                    min_proj = proj[i];
                    m2 = i;
                }
            swap(m2, m1 + 1);
            u = (proj[m1] + proj[m1 + 1]) / 2;
            MPI_Send(proj + m1 + 1, sum1 - (m1 + 1), MPI_DOUBLE, id + 1, 0, comm);
            MPI_Send(p_aux + (m1 + 1) * n_dims, (sum1 - (m1 + 1)) * n_dims, MPI_DOUBLE, id + 1, 1, comm);
            MPI_Send(idx_global + m1 + 1, sum1 - (m1 + 1), MPI_LONG, id + 1, 2, comm);
            sum1 = m1 + 1;
        } else {
            u = proj[m1];
            MPI_Send(proj + m1, sum1 - m1, MPI_DOUBLE, id + 1, 0, comm);
            MPI_Send(p_aux + m1 * n_dims, (sum1 - m1) * n_dims, MPI_DOUBLE, id + 1, 1, comm);
            MPI_Send(idx_global + m1, sum1 - m1, MPI_LONG, id + 1, 2, comm);
            sum1 = m1;
        }
    } else if ((id == p / 2 + 1) && (p % 2 == 1)) {
        u = 0;
        MPI_Recv(proj + sum1, total_size - (n_points / 2), MPI_DOUBLE, id - 1, 0, comm, MPI_STATUS_IGNORE);
        MPI_Recv(p_aux + sum1 * n_dims, (total_size - (n_points / 2)) * n_dims, MPI_DOUBLE, id - 1, 1, comm, MPI_STATUS_IGNORE);
        MPI_Recv(idx_global + sum1, total_size - (n_points / 2), MPI_LONG, id - 1, 2, comm, MPI_STATUS_IGNORE);
        sum1 += total_size - (n_points / 2);
    } else if ((n_points % 2 == 1) && (p % 2 == 0)) {
        if (id == p / 2 - 1) {
            if (total_size > n_points / 2) {
                if (total_size - sum1 >= n_points / 2) {
                    fprintf(stderr, "MEDIAN OUT OF BOUNDS - m1: %ld, sum1: %d\n", (n_points) / 2 - (total_size - sum1), sum1);
                    fflush(stderr);
                }
                m1 = quickselect_seq((n_points) / 2 - (total_size - sum1), 0, sum1);
                u = proj[m1];
                MPI_Send(proj + m1, sum1 - m1, MPI_DOUBLE, id + 1, 0, comm);
                MPI_Send(p_aux + m1 * n_dims, (sum1 - m1) * n_dims, MPI_DOUBLE, id + 1, 1, comm);
                MPI_Send(idx_global + m1, sum1 - m1, MPI_LONG, id + 1, 2, comm);
                sum1 = m1;
            } else {
                u = 0;
                MPI_Recv(proj + sum1, (n_points / 2) - total_size, MPI_DOUBLE, id + 1, 0, comm, MPI_STATUS_IGNORE);
                MPI_Recv(p_aux + sum1 * n_dims, ((n_points / 2) - total_size) * n_dims, MPI_DOUBLE, id + 1, 1, comm, MPI_STATUS_IGNORE);
                MPI_Recv(idx_global + sum1, (n_points / 2) - total_size, MPI_LONG, id + 1, 2, comm, MPI_STATUS_IGNORE);
                sum1 += (n_points / 2) - total_size;
            }
        } else if (id == p / 2) {
            if (total_size <= n_points / 2) {
                if (total_size + sum1 <= n_points / 2) {
                    fprintf(stderr, "MEDIAN OUT OF BOUNDS - m1: %ld, sum1: %d\n", (n_points) / 2 - (total_size), sum1);
                    fflush(stdout);
                }
                m1 = quickselect_seq((n_points) / 2 - (total_size), 0, sum1);
                u = proj[m1];
                MPI_Send(proj, m1, MPI_DOUBLE, id - 1, 0, comm);
                MPI_Send(p_aux, m1 * n_dims, MPI_DOUBLE, id - 1, 1, comm);
                MPI_Send(idx_global, m1, MPI_LONG, id - 1, 2, comm);
                sum1 = sum1 - m1;
                for (int i = 0; i < sum1; i++) {
                    proj[i] = proj[i + m1];
                    idx_global[i] = idx_global[i + m1];
                    for (int j = 0; j < n_dims; j++) {
                        p_aux[i * n_dims + j] = p_aux[(i + m1) * n_dims + j];
                    }
                }
            } else {
                u = 0;
                MPI_Recv(proj + sum1, total_size - (n_points / 2), MPI_DOUBLE, id - 1, 0, comm, MPI_STATUS_IGNORE);
                MPI_Recv(p_aux + sum1 * n_dims, (total_size - (n_points / 2)) * n_dims, MPI_DOUBLE, id - 1, 1, comm, MPI_STATUS_IGNORE);
                MPI_Recv(idx_global + sum1, total_size - (n_points / 2), MPI_LONG, id - 1, 2, comm, MPI_STATUS_IGNORE);
                sum1 += total_size - (n_points / 2);
            }
        } else {
            u = 0;
        }
    } else if ((n_points % 2 == 0) && (p % 2 == 0)) {
        if (id == p / 2 - 1) {
            if (total_size >= n_points / 2) {
                if (total_size - sum1 >= n_points / 2) {
                    fprintf(stderr, "MEDIAN OUT OF BOUNDS - m1: %ld, sum1: %d\n", (n_points) / 2 - 1 - (total_size - sum1), sum1);
                    fflush(stderr);
                }
                m1 = quickselect_seq((n_points) / 2 - 1 - (total_size - sum1), 0, sum1);
                if (total_size != n_points / 2) {
                    min_proj = proj[m1 + 1];
                    m2 = m1 + 1;
                    for (int i = m1 + 1; i < sum1; i++)
                        if (proj[i] < min_proj) {
                            min_proj = proj[i];
                            m2 = i;
                        }
                    swap(m2, m1 + 1);
                    u = (proj[m1] + proj[m1 + 1]) / 2;
                    MPI_Send(proj + m1 + 1, sum1 - (m1 + 1), MPI_DOUBLE, id + 1, 0, comm);
                    MPI_Send(p_aux + (m1 + 1) * n_dims, (sum1 - (m1 + 1)) * n_dims, MPI_DOUBLE, id + 1, 1, comm);
                    MPI_Send(idx_global + m1 + 1, sum1 - (m1 + 1), MPI_LONG, id + 1, 2, comm);
                    sum1 = m1 + 1;
                } else {
                    u = proj[m1] / 2;
                }
            } else {
                u = 0;
                MPI_Recv(proj + sum1, (n_points / 2) - total_size, MPI_DOUBLE, id + 1, 0, comm, MPI_STATUS_IGNORE);
                MPI_Recv(p_aux + sum1 * n_dims, ((n_points / 2) - total_size) * n_dims, MPI_DOUBLE, id + 1, 1, comm, MPI_STATUS_IGNORE);
                MPI_Recv(idx_global + sum1, (n_points / 2) - total_size, MPI_LONG, id + 1, 2, comm, MPI_STATUS_IGNORE);
                sum1 += (n_points / 2) - total_size;
            }
        } else if (id == p / 2) {
            if (total_size < n_points / 2) {
                if (n_points / 2 - 1 - total_size >= sum1 - 1) {
                    fprintf(stderr, "MEDIAN OUT OF BOUNDS - m1: %ld, sum1: %d\n", n_points / 2 - 1 - total_size, sum1);
                    fflush(stdout);
                }
                m1 = quickselect_seq(n_points / 2 - 1 - total_size, 0, sum1);
                min_proj = proj[m1 + 1];
                m2 = m1 + 1;
                // printf("HERE0 - m1: %d, m2: %d, proj[m1]: %f, proj[m2]: %f, min_proj: %f\n", m1, m2, proj[m1], proj[m2], min_proj);
                // fflush(stdout);
                for (int i = m1 + 2; i < sum1; i++) {
                    // printf("HERE0.1 - i: %d proj[i]: %f, min_proj: %f pt ", i, proj[i], min_proj);
                    // print_point(&p_aux[i*n_dims], n_dims);
                    // fflush(stdout);
                    if (proj[i] < min_proj) {
                        // printf("HERE1 - i: %d proj[i]: %f, pt ", i, proj[i]);
                        // print_point(&p_aux[i*n_dims], n_dims);
                        // fflush(stdout);
                        min_proj = proj[i];
                        m2 = i;
                    }
                }
                swap(m2, m1 + 1);
                // printf("HERE - m1: %d proj[m1]: %f, m1+1: %d, m2: %d proj[m2]: %f", m1, proj[m1], m1+1,m2,proj[m1+1]);
                // for (int j = 0; j < n_dims; j++) {
                //     printf(" %f ", p_aux[m1*n_dims+j]);
                //     fflush(stdout);
                // }
                // printf("\n");
                // for (int j = 0; j < n_dims; j++) {
                //     printf(" %f ", p_aux[(m1+1)*n_dims+j]);
                //     fflush(stdout);
                // }
                // printf("\n");
                u = (proj[m1] + proj[m1 + 1]) / 2;
                MPI_Send(proj, m1 + 1, MPI_DOUBLE, id - 1, 0, comm);
                MPI_Send(p_aux, (m1 + 1) * n_dims, MPI_DOUBLE, id - 1, 1, comm);
                MPI_Send(idx_global, m1 + 1, MPI_LONG, id - 1, 2, comm);
                sum1 = sum1 - (m1 + 1);
                for (int i = 0; i < sum1; i++) {
                    proj[i] = proj[i + m1 + 1];
                    idx_global[i] = idx_global[i + m1 + 1];
                    for (int j = 0; j < n_dims; j++) {
                        p_aux[i * n_dims + j] = p_aux[(i + m1 + 1) * n_dims + j];
                    }
                }
            } else if (total_size == n_points / 2) {
                m1 = -1;
                min_proj = proj[m1 + 1];
                m2 = m1 + 1;
                for (int i = m1 + 1; i < sum1; i++)
                    if (proj[i] < min_proj) {
                        min_proj = proj[i];
                        m2 = i;
                    }
                swap(m2, m1 + 1);
                u = proj[0] / 2;
            } else {
                u = 0;
                MPI_Recv(proj + sum1, total_size - (n_points / 2), MPI_DOUBLE, id - 1, 0, comm, MPI_STATUS_IGNORE);
                MPI_Recv(p_aux + sum1 * n_dims, (total_size - (n_points / 2)) * n_dims, MPI_DOUBLE, id - 1, 1, comm, MPI_STATUS_IGNORE);
                MPI_Recv(idx_global + sum1, total_size - (n_points / 2), MPI_LONG, id - 1, 2, comm, MPI_STATUS_IGNORE);
                sum1 += total_size - (n_points / 2);
            }
        }
    } else {
        u = 0;
    }

    block_size[id_initial] = sum1;

    if (sum1 > aux || sum1 <= 0) {
        fprintf(stderr, "ERROR! - id: %d, p: %d, sum: %d, aux: %d\n", id, p, sum1, aux);
        fflush(stdout);
        exit(0);
    }

    /* ALL PROCESSORS GET THE SCALAR PROJECTION OF THE MEDIAN POINT */
    double u_aux;
    MPI_Allreduce(&u, &u_aux, 1, MPI_DOUBLE, MPI_SUM, comm);
    u = u_aux;

    if (print_info == 1) {
        MPI_Barrier(comm);
        printf("id: %d, sum1: %d, u: %f\n", id, sum1, u);
        fflush(stdout);
        MPI_Barrier(comm);
    }

    /* ALL PROCESSORS COMPUTE THE CENTER POINT */
    double abnorm = 0.0;
    for (int i = 0; i < n_dims; ++i) {
        u_aux = (pt_arr_b[i] - pt_arr_a[i]);
        center[i] = u * u_aux;
        abnorm += u_aux * u_aux;
    }
    for (int i = 0; i < n_dims; ++i) {
        center[i] = center[i] / abnorm + pt_arr_a[i];
    }

    /* ALL PROCESSORS COMPUTE THE LOCAL RADIUS */
    double max_d = 0;
    for (int i = 0; i < block_size[id_initial]; i++) {
        d = 0;
        for (int j = 0; j < n_dims; j++)
            d += (center[j] - p_aux[i * n_dims + j]) * (center[j] - p_aux[i * n_dims + j]);
        if (d > max_d)
            max_d = d;
    }

    double rad = 0;

    /* THE GLOBAL RADIUS IS COMPUTED */
    MPI_Reduce(&max_d, &rad, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

    long left = -1;
    long right = -1;

    /* CHECKS IF THE PROCESSES ARE IN THE LAST LEVEL OF MPI PARALLELISM. IF SO, THE IDs ARE COMPUTED ACCORDINGLY */
    if (lvl == log(p_initial) / log(2) - 1) {
        /* EACH PROCESS LETS EACH OTHER KNOW OF ITS CURRENT SIZE */
        MPI_Allgather(&sum1, 1, MPI_INT, block_size, 1, MPI_INT, MPI_COMM_WORLD);
        // branch_size[0] = 2 * block_size[0] - 1;
        // for (int i = 1; i < p_initial; i++) {
        //     branch_size[i] = 2 * block_size[i] - 1 + branch_size[i - 1];
        // }
        // int first_id = pow(2, lvl + 1) - 1;
        // if (id_initial < 2) {
        //     left = first_id;
        //     right = first_id + branch_size[0];
        // } else if (id == 0) {
        //     left = first_id + branch_size[id_initial - 1];
        //     right = first_id + branch_size[id_initial];
        // } else {
        //     left = first_id + branch_size[id_initial - 2];
        //     right = first_id + branch_size[id_initial - 1];
        // }
    // } else {
    }
        left = 2 * n_id + 1;
        right = 2 * n_id + 2;
    // }

    if (print_info == 1) {
        printf("id: %d, id_initial: %d, left: %ld, right: %ld\n", id, id_initial, left, right);
        fflush(stdout);
    }

    /* ONLY THE PROCESS WITH THE CURRENT RANK 0 SAVES THE TREE NODE CORRESPONDING TO THIS ballAlg CALL */
    if (!id) {
        tree[tree_counter].node_id = n_id;
        tree[tree_counter].center_idx = n_center;
        tree[tree_counter].radius = sqrt(rad);
        tree[tree_counter].left = left;
        tree[tree_counter].right = right;
        for (int i = 0; i < n_dims; ++i) {
            centers[n_center][i] = center[i];
        }
        n_center++;
        tree_counter++;
    }

    /* LAST MPI PARELLEL LEVEL */
    if (lvl == log(p_initial) / log(2) - 1) {

        free(pt_arr_a);
        free(pt_arr_b);

        /* STRUCTURE FOR THE SEQUENTIAL VERSION */
        pts = (pt *)malloc(sizeof(pt) * block_size[id_initial]);
        for (int i = 0; i < block_size[id_initial]; ++i) {
            pts[i].pt = &p_aux[i * n_dims];
            pts[i].i = idx_global[i];
        }
        /* ONE PROCESS (LOWER RANK) COMPUTES THE LEFT BRANCH... */
        if (threads == 1 || (threads & (threads - 1)) != 0)
            max_parallel_level = 0;
        else
            max_parallel_level = ceil(log(threads) / log(2)) - 1;
        id_last = pow(2, ceil(log(np)/log(2))) - 1; // First index of the last level
        int sum = 0;
        n_levels = ceil(log(block_size[0]) / log(2)); // Number of levels in the tree minus one
        for (int i = 1; i <= id_initial; i++) {
            n_levels = ceil(log(block_size[i-1]) / log(2)); // Number of levels in the tree minus one
            sum += 2 * block_size[i-1] - pow(2, n_levels);
            id_last += 2 * block_size[i-1] - pow(2, n_levels);
        }

        sum+=2 * block_size[p_initial-1] - pow(2, n_levels);
        // printf("id: %d, sum: %d, exp: %d\n", id_initial, sum, (int)(2*np-pow(2, ceil(log(np)/log(2)))));
        // fflush(stdout);

        if (!id) {
            // tree_idx = left;
            // printf("id: %d, id_last: %ld\n", id_initial, id_last);
            ballAlg_omp(0, block_size[id_initial], left, 0, threads);
            /* ... AND THE OTHER THE RIGHT BRANCH */
        } else {
            // tree_idx = right;
            // printf("id: %d, id_last: %ld\n", id_initial, id_last);
            ballAlg_omp(0, block_size[id_initial], right, 0, threads);
        }
        /* NOT YET IN THE LAST MPI-PARALLEL LEVEL */
    } else {
        /* SPLITS THE PROCESSES IN TWO: HALF WILL BE DOING THE LEFT BRANCH, THE REMAINING THE RIGHT BRANCH */
        MPI_Comm new_comm;
        MPI_Comm_split(comm, (id + p % 2) / (p / 2 + p % 2), id, &new_comm);

        if (n_id != 0)
            MPI_Comm_free(&comm);

        MPI_Comm_rank(new_comm, &new_id);
        MPI_Comm_size(new_comm, &p);

        if (print_info == 1) {
            printf("id: %d, new_id: %d, n_id: %ld, proc: %d\n", id, new_id, n_id, p);
        }

        if (id < p) {
            id = new_id;
            ballAlg_mpi(n_points / 2, left, lvl + 1, new_comm, threads);
        } else {
            id = new_id;
            ballAlg_mpi(n_points / 2 + n_points % 2, right, lvl + 1, new_comm, threads);
        }
    }
}

void print_tree(node *tree) {

    if (!id_initial) {
        fprintf(stdout, "%d %ld\n", n_dims, 2 * np - 1);
        fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (long i = 0; i < tree_counter; ++i) {
        fprintf(stdout, "%ld %ld %ld %f ", tree[i].node_id, tree[i].left, tree[i].right, tree[i].radius);
        // fprintf(stdout, "%f ", tree[i].radius);
        fflush(stdout);
        if (tree[i].left == -1) {
            print_point(pts[tree[i].center_idx].pt, n_dims);
            fflush(stdout);
        } else {
            print_point(centers[tree[i].center_idx], n_dims);
            fflush(stdout);
        }
    }
}

int main(int argc, char **argv) {

    double elapsed_time;
    int n_threads;

    MPI_Init(&argc, &argv);

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    omp_set_nested(1);

    if (argc != 4) {
        if (!id) {
            printf("Usage: %s <n_dims> <np> <seed>\n", argv[0]);
            fflush(stdout);
        }
        MPI_Finalize();
        exit(1);
    }

    n_dims = atoi(argv[1]);
    np = atol(argv[2]);
    unsigned seed = atoi(argv[3]);
    srandom(seed);

    p_initial = p;
    id_initial = id;

    aux = ceil(np / p) * 5 / 2;
    p_aux = (double *)malloc(n_dims * aux * sizeof(double));
    // pt_arr = (double **)malloc(aux * sizeof(double *));
    // for (long i = 0; i < aux; i++)
    //     pt_arr[i] = &p_aux[i * n_dims];

    block_size = (int *)malloc(p * sizeof(int));
    branch_size = (int *)malloc(p * sizeof(int));
    for (int i = 0; i < p; i++)
        block_size[i] = BLOCK_SIZE(i, p, np);

    idx_global = (long *)malloc(aux * sizeof(long));

    int BL = BLOCK_LOW(id, p, np);
    int BH = BLOCK_HIGH(id, p, np);
    for (int i = 0; i < np; i++)
        for (int j = 0; j < n_dims; j++) {
            if ((i >= BL) && (i <= BH)) {
                p_aux[(i - BL) * n_dims + j] = RANGE * ((double)random()) / RAND_MAX;
            } else
                random();
        }

    pt_arr_idx0 = (double *)malloc(n_dims * sizeof(double));
    pt_arr_a = (double *)malloc(n_dims * sizeof(double));
    pt_arr_b = (double *)malloc(n_dims * sizeof(double));

    proj = (double *)malloc(aux * sizeof(double));

    pivots = (double *)malloc(p * sizeof(double));
    glob_pivots = (double *)malloc(p * p * sizeof(double));
    send_displs = (int *)malloc(p * sizeof(int));
    recv_displs = (int *)malloc(p * sizeof(int));
    send_counts = (int *)malloc(p * sizeof(int));
    recv_counts = (int *)malloc(p * sizeof(int));

    recv_buffer = (double *)malloc(aux * sizeof(double));
    recv_buffer_long = (long *)malloc(aux * sizeof(long));
    p_aux_2 = (double *)malloc(n_dims * aux * sizeof(double));

    sizes = (int *)malloc(p * sizeof(int));

    tree = (node *)malloc((ceil(log(p_initial) / log(2)) + (2 * aux - 1)) * sizeof(node));
    double *center_aux = (double *)malloc((ceil(log(p_initial) / log(2)) + aux - 1) * n_dims * sizeof(double));
    centers = (double **)malloc((ceil(log(p_initial) / log(2)) + aux - 1) * sizeof(double *));
    center = (double *)malloc(n_dims * sizeof(double));

#pragma omp parallel
    {
#pragma omp single
        n_threads = omp_get_num_threads();
#pragma omp for
        for (int i = 0; i < block_size[id]; i++)
            idx_global[i] = BLOCK_LOW(id, p, np) + i;
        int max_it = ceil(log(p_initial) / log(2)) + aux - 1;
#pragma omp for
        for (int i = 0; i < max_it; ++i)
            centers[i] = &center_aux[i * n_dims];
    }

    ballAlg_mpi(np, 0, 0, MPI_COMM_WORLD, n_threads);

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();

    if (!id_initial) {
        fprintf(stderr, "%.1lf\n", elapsed_time);
        fflush(stdout);
    }

    print_tree(tree);

    free(p_aux);
    free(block_size);
    free(branch_size);
    free(idx_global);
    free(pt_arr_idx0);
    free(proj);
    free(pivots);
    free(glob_pivots);
    free(send_displs);
    free(recv_displs);
    free(send_counts);
    free(recv_counts);
    free(recv_buffer);
    free(recv_buffer_long);
    free(p_aux_2);
    free(tree);
    free(center_aux);
    free(centers);
    free(center);
    free(pts);
    free(sizes);

    MPI_Finalize();
}