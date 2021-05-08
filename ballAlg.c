#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "gen_points.h"

#define BLOCK_LOW(id,p,np) ((id)*(np)/(p))
#define BLOCK_HIGH(id,p,np) (BLOCK_LOW((id)+1,p,np)-1)
#define BLOCK_SIZE(id,p,np) (BLOCK_HIGH(id,p,np)-BLOCK_LOW(id,p,np)+1)

typedef struct tree_node {
    int center_idx;
    double radius;
    long left;
    long right;
} node;

typedef struct pt_struct {
    int i;
    double proj;
    double *pt;
} pt;

pt *pts;

double **centers;
long n_center = 0;

int n_dims; // Dimensions of the input points
long np;    // Number of generated points
/* Array of points: allocated at the beggining of the program;
* Points are randomly generated at the start of the program and don't change as it runs
*/
double **pt_array;
/* Struct that saves all tree information, for later printing;
* Has 2 * np - 1 points
*/
node *tree;
/* vector with the sorted indices. 
* All the computations are made in this vector:
* sorting according to projection value, splitting in right/left ball
* ...
*/

/* Function to print a point to a given file */
void print_point(double *pt, int dim, FILE *fd) {
    for (int i = 0; i < dim; ++i)
        fprintf(fd, "%f ", pt[i]);
    fprintf(fd, "\n");
}

/* Computes the squared distance between 2 points */
double dist(double *pt1, double *pt2) {
    double dist = 0.0;
    for (int d = 0; d < n_dims; ++d)
        dist += (pt1[d] - pt2[d]) * (pt1[d] - pt2[d]);
    return dist;
}

void swap(long *i, long *j) {
    long temp = *i;
    *i = *j;
    *j = temp;
}

void swap_pt(long i, long j) {
    pt temp = pts[i];
    pts[i] = pts[j];
    pts[j] = temp;
}

void furthest_apart(int l, int r, long *a_idx, long *b_idx) {
    double max_d = 0.0;
    double d;
    long idx0 = pts[l].i;
    for (int i = l + 1; i < r; ++i) {
        if (pts[l].i < idx0)
            idx0 = l;
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
        swap(a_idx, b_idx);
    *a_idx = pts[*a_idx].i;
    *b_idx = pts[*b_idx].i;
}

double orth_projv1(double *a, double *b, double *p) {
    // Orthogonal projection in the first coordinate
    double proj = 0.0;
    for (int i = 0; i < n_dims; ++i)
        proj += (p[i] - a[i]) * (b[i] - a[i]);

    return proj;
}

void orth_projv2(double *a, double *b, long idx_p, double *ab_proj) {
    double abnorm = 0.0, aux, u = pts[idx_p].proj;
    for (int i = 0; i < n_dims; ++i) {
        aux = (b[i] - a[i]);
        ab_proj[i] = u * aux;
        abnorm += aux * aux;
    }
    for (int i = 0; i < n_dims; ++i)
        ab_proj[i] = ab_proj[i] / abnorm + a[i];
}

int is_seq(int i, int j) {
    return pts[i].proj <= pts[j].proj;
}

/* Gets the k-th element of a piece of the projections array (from l to h):
* changes only the idx (indices) array;
* the comparison critirion is the first coordinate of all projections (stored in proj_scalar);
* also partitions the data: at the left of the k-th element of the indices array are all the indices
*   corresponding to points whose the first coordinate projection in the current ab line is not bigger
*   than that of the k-th element (points that should be on the left ball);
* at the right of the k-th element of the indices array are all the indices
*   corresponding to points whose the first coordinate projection in the current ab line is not smaller
*   than that of the k-th element (points that should be on the right ball).
*/
int get_kth_element(int k, int l, int h) {

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
        *median_2 = *median_1;
        for (int i = l + (h - l) / 2; i < h; ++i) {
            if (is_seq(*median_2, i)) {
                *median_2 = i;
            }
        }
        // *median_2 = get_kth_element(0, l + (h - l) / 2, h);
    }
}

double *abproj;

long n_nodes; // Total number of tree nodes (in the end will equal 2 * np - 1)
long n_levels;
long id_last;

int id, p;
void ballAlg(long l, long r, long tree_id, int lvl) {

    for (int i = 0; i < BLOCK_SIZE(id,p,np); ++i) {
        printf("id: %d, pontos:", id);
        print_point(pts[i].pt, n_dims, stdout);
        printf("id: %d, i: %d\n", id, pts[i].i);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    exit(0);

    // printf("l: %ld, r: %ld, tree_id: %ld, lvl: %d\n", l, r, tree_id, lvl);

    if (r - l == 1) {
        tree[tree_id].center_idx = l;
        tree[tree_id].radius = 0;
        tree[tree_id].left = -1;
        tree[tree_id].right = -1;
        return;
    }

    long c_id = n_center++;

    long a, b;
    // 2. Compute points a and b, furthest apart in the current set (approx)
    furthest_apart(l, r, &a, &b);
    printf("a: %ld, b: %ld\n", a, b);

    // 3. Perform the orthogonal projection of all points onto line ab
    for (int i = l; i < r; ++i)
        pts[i].proj = orth_projv1(pt_array[a], pt_array[b], pts[i].pt);

    // 4. Compute the center, defined as the median point over all projections
    int m1, m2 = -1;
    get_median(l, r, &m1, &m2);
    printf("m1: %d, m2: %d\n", m1, m2);

    if ((r - l) % 2) {
        orth_projv2(pt_array[a], pt_array[b], m1, centers[c_id]);
    } else {
        double abnorm = 0.0, aux, u = (pts[m1].proj + pts[m2].proj) / 2;
        for (int i = 0; i < n_dims; ++i) {
            aux = (pt_array[b][i] - pt_array[a][i]);
            centers[c_id][i] = u * aux;
            abnorm += aux * aux;
        }
        for (int i = 0; i < n_dims; ++i)
            centers[c_id][i] = centers[c_id][i] / abnorm + pt_array[a][i];
    }

    double max_r = 0, rad;
    for (int i = l; i < r; ++i) {
        rad = dist(centers[c_id], pts[i].pt);
        if (rad > max_r)
            max_r = rad;
    }
    tree[tree_id].center_idx = c_id;
    tree[tree_id].radius = sqrt(max_r);

    if (((np & (np - 1)) != 0) && lvl == n_levels - 1) {
        tree[tree_id].left = id_last;
        tree[tree_id].right = id_last + 1;
        id_last += 2;
    } else {
        tree[tree_id].left = 2 * tree_id + 1;
        tree[tree_id].right = 2 * tree_id + 2;
    }

    ballAlg(l, l + (r - l) / 2, tree[tree_id].left, lvl + 1);
    ballAlg(l + (r - l) / 2, r, tree[tree_id].right, lvl + 1);
}

void print_tree(node *tree) {
    fprintf(stdout, "%d %ld\n", n_dims, 2 * np - 1);
    for (long i = 0; i < 2 * np - 1; ++i) {
        fprintf(stdout, "%ld %ld %ld %f ", i, tree[i].left, tree[i].right, tree[i].radius);
        if (tree[i].left == -1)
            print_point(pts[tree[i].center_idx].pt, n_dims, stdout);
        else
            print_point(centers[tree[i].center_idx], n_dims, stdout);
    }
}

int main(int argc, char **argv) {

    double elapsed_time;

    MPI_Init (&argc, &argv);

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);

    if (argc != 4) {
        if (!id)
            printf ("Command line: %s <m>\n", argv[0]);
        MPI_Finalize();
        exit (1);
    }

    pt_array = get_points(argc, argv, &n_dims, &np);

    pts = (pt *)malloc(sizeof(pt) * ceil(np/p));
    
    n_nodes = 2 * np - 1;
    n_levels = ceil(log(np) / log(2));
    id_last = pow(2, n_levels) - 1;

    for (int i = 0; i < BLOCK_SIZE(id,p,np); ++i) {
        pts[i].pt = pt_array[BLOCK_LOW(id,p,np)+i];
        pts[i].i = BLOCK_LOW(id,p,np)+i;
    }

    tree = (node *)malloc((2 * np - 1) * sizeof(node));
    double *center_aux = (double *)malloc((np - 1) * n_dims * sizeof(double));
    centers = (double **)malloc((np - 1) * sizeof(double *));
    for (int i = 0; i < np - 1; ++i)
        centers[i] = &center_aux[i * n_dims];

    elapsed_time += MPI_Wtime();
    if (!id) {
        printf ("Total elapsed time: %10.6f\n", elapsed_time);
    }

    ballAlg(0, np, 0, 0);

    MPI_Finalize ();

    // fprintf(stderr, "%.1lf\n", elapsed_time);

    print_tree(tree);

    free(pt_array[0]);
    free(pt_array);

    free(centers[0]);
    free(centers);

    free(tree);

    free(pts);

    exit(0);
}
