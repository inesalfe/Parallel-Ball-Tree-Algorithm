#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "gen_points.h"
// #define DEBUG

typedef struct tree_node {
    double *center;
    double radius;
    long left;
    long right;
} node;

int n_dims; // Dimensions of the input points
long np;    // Number of generated points
/* Array of points: allocated at the beggining of the program;
* Points are randomly generated at the start of the program and don't change as it runs.
*/
double **pt_array;
/* Struct that saves all tree information, for later printing;
* Has 2 * np - 1 points.
*/
node *tree;
/* Array with the orthogonal projections onto line ab, at a given time;
* Only one array is needed since the projections from seperate sets are non-overlapping
* (that is, at any given time, each point will only be projected in a unique line ab, 
* even if a and b differ from subset to subset). */
double *proj_scalar;
/* Array with the sorted indices. 
* All the computations are made in this vector:
* sorting according to projection value, splitting in right/left ball
* ...
*/
long *idx;

#pragma omp threadprivate(n_dims, np, pt_array, tree, proj_scalar)

/* Function to print a point to a given file */
void print_point(double *pt, int dim, FILE *fd) {
    for (int i = 0; i < dim; ++i)
        fprintf(fd, "%f ", pt[i]);
    fprintf(fd, "\n");
}

/* Computes the squared distance between 2 points;
* since the ordering of the values is the same with or without the square root,
* and most times we are only interested in the relative order and not the actual 
* distance value, the square root isn't taken, in order to improve efficiency.
* Everytime the actual distance is needed (when computing the radius of a set), 
* the square root is computed outside the scope of this functions. */
double dist(double *pt1, double *pt2) {
    double dist = 0.0;
    for (int d = 0; d < n_dims; ++d)
        dist += (pt1[d] - pt2[d]) * (pt1[d] - pt2[d]);
    return dist;
}

/* Function to swap two values */
void swap(long *i, long *j) {
    long temp = *i;
    *i = *j;
    *j = temp;
}

/* Computes an approximation of the two furthest apart points from the set */
void furthest_apart(int l, int r, long *a_idx, long *b_idx) {
    double max_d = 0.0;
    double d;
    /* Finds the element with lower index in the current set, which might not be
    * idx[l];
    * Not necessarily needed for a correct algorithm, but helps reproduce
    * the teaching staff's tree, for easier comparison. */
    long idx0 = idx[l];
    for (int i = l + 1; i < r; ++i) {
        if (idx[i] < idx0)
            idx0 = idx[i];
    }
    for (int i = l; i < r; ++i) {
        d = dist(pt_array[idx0], pt_array[idx[i]]);
        if (d > max_d) {
            max_d = d;
            *a_idx = idx[i];
        }
    }
    max_d = 0;
    for (int i = l; i < r; ++i) {
        d = dist(pt_array[*a_idx], pt_array[idx[i]]);
        if (d > max_d) {
            max_d = d;
            *b_idx = idx[i];
        }
    }
    /* Point a has always a smaller first coordinate than b (again, for easier comparison
    * of results with those of the teaching staff). */
    if (pt_array[*a_idx][0] > pt_array[*b_idx][0])
        swap(a_idx, b_idx);
}

/* Orthogonal projection onto line ab - scalar:
* computes the unnormalized length of the projected segment, with respect to point a.
* May wield negative values if a point is at the left of a, which may occur since ab is an
* approximation. This isn't problematic, since the ordering will still be valid.
* The formula used is (p-a).(b-a), where p is the point from which we want to compute
* the projection, and . denotes the dot product. */
double orth_projv1(double *a, double *b, double *p) {
    // Orthogonal projection in the first coordinate
    double proj = 0.0;
    for (int i = 0; i < n_dims; ++i)
        proj += (p[i] - a[i]) * (b[i] - a[i]);

    return proj;
}

/* Orthogonal projection onto line ab - vector:
* Computes the actual orthogonal projection in all coordinates, and stores the resulting
* vector into ab_proj. Uses the value already computed by the previous function (orth_projv1). */
void orth_projv2(double *a, double *b, long idx_p, double *ab_proj) {
    double abnorm = 0.0, aux, u = proj_scalar[idx_p];
    // #pragma omp parallel for private(aux) reduction(+:abnorm)
        for (int i = 0; i < n_dims; ++i) {
            aux = (b[i] - a[i]);
            ab_proj[i] = u * aux;
            abnorm += aux * aux;
        }

// Como é que eu aproveito as threads criadas em cima para o for de baixo?
    // #pragma omp parallel for
        for (int i = 0; i < n_dims; ++i)
            ab_proj[i] = ab_proj[i] / abnorm + a[i];
}

/* Comparison function: compares the orthogonal projections at indexes i and j. */
int is_seq(int i, int j) {
    return proj_scalar[i] <= proj_scalar[j];
}

int cmp(const void *a, const void *b) {
    return (proj_scalar[*(long *)a] > proj_scalar[*(long *)b] ? 1 : -1);
}

/* Gets the k-th element of a piece of the projections array (from l to h):
* changes only the idx (indices) array;
* the comparison critirion is the scalar projection of all points in the set (stored in proj_scalar);
* also partitions the data: at the left of the k-th element of the indices array are all the indices
*   corresponding to points whose projection in the current ab line is not bigger
*   than that of the k-th element (points that should be on the left ball);
* at the right of the k-th element of the indices array are all the indices
*   corresponding to points whose the first coordinate projection in the current ab line is not smaller
*   than that of the k-th element (points that should be on the right ball).
*/
int get_kth_element(int k, int l, int h) {

    if (h == l + 1)
        return idx[l];

    else {
        int i = l;
        int j = h;

        while (i < j) {
            do {
                i++;
            } while (i < h && is_seq(idx[i], idx[l]));
            if (i == h) {
                j = h - 1;
                break;
            }
            do {
                j--;
            } while (!is_seq(idx[j], idx[l]));
            if (i < j)
                swap(&idx[i], &idx[j]);
        }

        swap(&idx[l], &idx[j]);

        if (j - l == k)
            return idx[j];
        if (k < j - l)
            return get_kth_element(k, l, j);
        else
            return get_kth_element(k - (j - l) - 1, j + 1, h);
    }
}

/* Same as before, but with a different partition method. Used for comparison. */
int get_kth_element2(int k, int l, int h) {

    if (h == l + 1)
        return idx[l];

    int i = l;
    long pivot = idx[l];

    for (int j = l + 1; j < h; ++j) {
        if (is_seq(idx[j], pivot)) {
            ++i;
            swap(&idx[i], &idx[j]);
        }
    }
    swap(&idx[i], &idx[l]);

    if (i - l == k)
        return idx[i];
    if (k < i - l)
        return get_kth_element2(k, l, i);
    else
        return get_kth_element2(k - (i - l) - 1, i + 1, h);
}

/* Gets the median point of the current set in N.log(N) (average), where N is the
* size of the set; Uses the partition method of quicksort, but only working with one of the
* sets. This way, it is always better, or, in the worst case, equal, than quicksort. */
void get_median(int l, int h, int *median_1, int *median_2) {
    if ((h - l) % 2 == 1)
        *median_1 = get_kth_element((h - l - 1) / 2, l, h);
    else {
        *median_1 = get_kth_element((h - l) / 2 - 1, l, h);
        *median_2 = get_kth_element(0, l + (h - l) / 2, h);
    }
}

/* Auxiliar arrays:
*  center1: when the set has an even number of points, stores one of the median points;
*  abproj: vector projection onto line ab. */
double *center1;
double *abproj;

long n_nodes = 0; // Total number of tree nodes (in the end will equal 2 * np - 1)

#pragma omp threadprivate(center1, abproj, n_nodes)

/* Actual algorithm to compute the tree. */
long ballAlg(long l, long r) {

    ++n_nodes;
    long id = n_nodes - 1;

#ifdef DEBUG
    printf("%ld %ld\n", l, r);
#endif

    if (r - l == 1) {
        for (int i = 0; i < n_dims; ++i)
            tree[id].center[i] = pt_array[idx[l]][i];
        tree[id].radius = 0;
        tree[id].left = -1;
        tree[id].right = -1;
        return id;
    }

    long a, b;
    // 2. Compute points a and b, furthest apart in the current set (approx)
    furthest_apart(l, r, &a, &b);

#ifdef DEBUG
    printf("a: %ld, b: %ld\n", a, b);
#endif

    // 3. Perform the orthogonal projection of all points onto line ab

    // #pragma omp parallel for
        for (int i = l; i < r; ++i)
            proj_scalar[idx[i]] = orth_projv1(pt_array[a], pt_array[b], pt_array[idx[i]]);

    // 4. Compute the center, defined as the median point over all projections
    int m1, m2 = -1;
    get_median(l, r, &m1, &m2);
    // qsort(idx + l, r - l, sizeof(long), cmp);

#ifdef DEBUG
    printf("Median index 1: %d, Median index 2: %d\nMedian point 1: ", m1, m2);
    print_point(pt_array[m1], n_dims, stdout);
#endif

#ifdef DEBUG
    printf("\tidx: ");
    for (int i = 0; i < np; ++i)
        printf("%ld ", idx[i]);
    printf("\n");
#endif

    if ((r - l) % 2) {
        orth_projv2(pt_array[a], pt_array[b], m1, tree[id].center);
    } else {
        orth_projv2(pt_array[a], pt_array[b], m1, tree[id].center);
        orth_projv2(pt_array[a], pt_array[b], m2, center1);
        for (int i = 0; i < n_dims; ++i)
            tree[id].center[i] = .5 * (tree[id].center[i] + center1[i]);
    }

    double max_r = 0, rad;
    for (int i = l; i < r; ++i) {
        rad = dist(tree[id].center, pt_array[idx[i]]);
        if (rad > max_r)
            max_r = rad;
    }
    max_r = sqrt(max_r);

    long L, R;

    L = ballAlg(l, l + (r - l) / 2);
    R = ballAlg(l + (r - l) / 2, r);

    tree[id].radius = max_r;
    tree[id].left = L;
    tree[id].right = R;

    return id;
}

/* Prints the resulting tree to a file or stdout. */
void print_tree(node *tree) {
    FILE *fd = fopen("pts.txt", "w");
    fprintf(fd, "%d %ld\n", n_dims, 2 * np - 1);
    for (long i = 0; i < 2 * np - 1; ++i) {
        fprintf(fd, "%ld %ld %ld %f ", i, tree[i].left, tree[i].right, tree[i].radius);
        print_point(tree[i].center, n_dims, fd);
    }
    fclose(fd);
}

/* Main: gets the initial points, allocates all memory and calls the routine that creates the tree (ballAlg);
* Also computes the time it takes for the above computations, printing it and the resulting tree afterwards;
* Lastly, frees all memory used. */
int main(int argc, char **argv) {

    double exec_time;
    exec_time = -omp_get_wtime();

    #pragma omp parallel num_threads(4)
    {

        #pragma omp single
        {
            // Isto aqui me baixo não funciona como devia...
            int total_n_threads = omp_get_num_threads();
            printf("Total number of threads: %d\n", total_n_threads);

            // 1. Get input sample points
            pt_array = get_points(argc, argv, &n_dims, &np);

            #ifdef DEBUG
                FILE *fd = fopen("points.txt", "w");
                for (int i = 0; i < np; ++i) {
                    print_point(pt_array[i], n_dims, stdout);
                    print_point(pt_array[i], n_dims, fd);
                }
                fclose(fd);
            #endif

            idx = (long *)malloc(sizeof(long) * np);

        }
    
        #pragma omp for
            for (int i = 0; i < np; ++i) {
                if (i == 0)
                    printf("Aqui!\n");
                idx[i] = i;
            }

        #pragma omp single
        {

            proj_scalar = (double *)malloc(np * sizeof(double));
            center1 = (double *)malloc(n_dims * sizeof(double));

            tree = (node *)malloc((2 * np - 1) * sizeof(node));
            for (int i = 0; i < (2 * np - 1); ++i)
                tree[i].center = (double *)malloc(n_dims * sizeof(double));

            ballAlg(0, np);

            exec_time += omp_get_wtime();
            fprintf(stderr, "%.3lf\n", exec_time);

            print_tree(tree);

            free(center1);

            for (int i = 0; i < (2 * np - 1); ++i)
                free(tree[i].center);
            free(tree);

            free(proj_scalar);
            free(idx);
            free(pt_array[0]);
            free(pt_array);

        }

    }

    return 0;
}
