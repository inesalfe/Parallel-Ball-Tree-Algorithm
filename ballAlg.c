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

// #pragma omp threadprivate(proj_scalar)

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
		for (int i = 0; i < n_dims; ++i) {
			aux = (b[i] - a[i]);
			ab_proj[i] = u * aux;
			abnorm += aux * aux;
		}
		for (int i = 0; i < n_dims; ++i)
		{
			ab_proj[i] = ab_proj[i] / abnorm + a[i];
		}
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
void get_kth_element2(int k, int l, int h) {

	if (h == l + 1)
		return;

	int i = l;
	long pivot = idx[l];

	for (int j = l + 1; j < h; ++j) {
		// for (int in = l; in < h; in++) {
		// 	printf("%ld ", idx[in]);
		// }
		// printf("\n");
		// for (int in = l; in < h; in++) {
		// 	printf("%f ", proj_scalar[idx[in]]);
		// }
		// printf("\n");
		if (is_seq(idx[j], pivot)) {
			++i;
			swap(&idx[i], &idx[j]);
		}		
	}
	
	swap(&idx[i], &idx[l]);

	if (i - l == k)
		return;

	if (k < i - l)
		// #pragma omp task
		get_kth_element2(k, l, i);

	else
		// #pragma omp task
		get_kth_element2(k - (i - l) - 1, i + 1, h);
}

/* Gets the median point of the current set in N.log(N) (average), where N is the
* size of the set; Uses the partition method of quicksort, but only working with one of the
* sets. This way, it is always better, or, in the worst case, equal, than quicksort. */
// void get_median(int l, int h, int *median_1, int *median_2) {

// 	if ((h - l) % 2 == 1)
// 		*median_1 = get_kth_element2((h - l - 1) / 2, l, h);
// 	else {
// 		*median_1 = get_kth_element2((h - l) / 2 - 1, l, h);
// 		*median_2 = get_kth_element2(0, l + (h - l) / 2, h);
// 	}
// }

/* Auxiliar arrays:
*  center1: when the set has an even number of points, stores one of the median points;
*  abproj: vector projection onto line ab. */
double *center1;
double *abproj;

long n_nodes; // Total number of tree nodes (in the end will equal 2 * np - 1)
long n_levels;
long index_last;
long n_last_level;
long l_lower;
long l_upper;

/* Actual algorithm to compute the tree. */
void ballAlg(long l, long r, long id) {

	printf("tid: %d\n", omp_get_thread_num());

	long a = -1;
	long b = -1;
	int m1, m2 = -1;
	double abnorm = 0.0;

	long idx0 = np;
	double max_d = 0.0;
	double max_r = 0.0;

	#pragma omp parallel num_threads(4)
	{

		if (r - l == 1) {
			#pragma omp for
			for (int i = 0; i < n_dims; ++i)
				tree[id].center[i] = pt_array[idx[l]][i];
			#pragma omp single
			{
				tree[id].radius = 0;
				tree[id].left = -1;
				tree[id].right = -1;
			}
		}

		else {

			long idx0_t = idx[l];
			#pragma omp for
				for (int i = l + 1; i < r; ++i) {
					if (idx[i] < idx0_t)
						idx0_t = idx[i];
				}

			#pragma omp critical
			{
				if (idx0_t < idx0)
					idx0 = idx0_t;
			}

			#pragma omp barrier

			double d;
			double max_d_t = 0.0;
			long a_t = -1;
			#pragma omp for
				for (int i = l; i < r; ++i) {
					d = dist(pt_array[idx0], pt_array[idx[i]]);
					if (d > max_d_t) {
						max_d_t = d;
						a_t = idx[i]; 
					}
				}

			#pragma omp critical
			{
				if (max_d_t > max_d) {
					max_d = max_d_t;
					a = a_t;
				}
			}

			#pragma omp barrier
			max_d = 0;

			long b_t = -1;
			max_d_t = 0;
			#pragma omp for
				for (int i = l; i < r; ++i) {
					d = dist(pt_array[a], pt_array[idx[i]]);
					if (d > max_d_t) {
							max_d_t = d;
							b_t = idx[i];
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
				if (pt_array[a][0] > pt_array[b][0])
					swap(&a, &b);
				// printf("%ld\n", a);
				// printf("%ld\n", b);
			}

			// 3. Perform the orthogonal projection of all points onto line ab
			#pragma omp for
				for (int i = l; i < r; ++i) {
					proj_scalar[idx[i]] = orth_projv1(pt_array[a], pt_array[b], pt_array[idx[i]]);
				}

			// #pragma omp single
			// {
			// 	for (int i = l; i < r; ++i) {
			// 		printf("idx: %ld, proj: %f\n", idx[i], proj_scalar[idx[i]]);
			// 	}
			// 	printf("\n");
			// }

			// 4. Compute the center, defined as the median point over all projections
			#pragma omp single
			{
				if ((r - l) % 2 == 1) {
					get_kth_element2((r - l - 1) / 2, l, r);
					m1 = idx[l + (r - l - 1) / 2];
				}
				else {
					get_kth_element2((r - l) / 2 - 1, l, r);
					get_kth_element2(0, l + (r - l) / 2, r);
					m1 = idx[l + (r - l) / 2 - 1];
					m2 = idx[l + (r - l) / 2];
				}
			}

			// #pragma omp single
			// {
			// 	for (int i = l; i < r; ++i) {
			// 		printf("idx: %ld, proj: %f\n", idx[i], proj_scalar[idx[i]]);
			// 	}
			// 	printf("m1: %d, m2: %d\n", m1, m2);
			// 	printf("\n");
			// }

			double aux, u = proj_scalar[m1];
			#pragma omp for reduction(+:abnorm)
				for (int i = 0; i < n_dims; ++i) {
					aux = (pt_array[b][i] - pt_array[a][i]);
					tree[id].center[i] = u * aux;
					abnorm += aux * aux;
				}

			if ((r - l) % 2 == 0) {
				u = proj_scalar[m2];
				#pragma omp for
					for (int i = 0; i < n_dims; ++i) {
						aux = (pt_array[b][i] - pt_array[a][i]);
						center1[i] = u * aux / abnorm + pt_array[a][i];
						tree[id].center[i] = 0.5 * (tree[id].center[i] / abnorm + pt_array[a][i] + center1[i]);
					}
			}

			else {
				#pragma omp for
					for (int i = 0; i < n_dims; ++i)
						tree[id].center[i] = tree[id].center[i] / abnorm + pt_array[a][i];
			}

			double rad;
			double max_r_t = 0;

			#pragma omp for
			for (int i = l; i < r; ++i) {
				rad = dist(tree[id].center, pt_array[idx[i]]);
				if (rad > max_r_t)
					max_r_t = rad;
			}

			#pragma omp critical
			{
				if (max_r_t > max_r)
					max_r = max_r_t;
			}
			#pragma omp barrier

			#pragma omp single
			{
				tree[id].radius = sqrt(max_r);
				if (((np & (np - 1)) != 0) && id >= l_lower && id <= l_upper) {
					tree[id].left = id + n_last_level - 1;
					tree[id].right = id + n_last_level;
				}
				else {
					tree[id].left = 2*id+1;
					tree[id].right = 2*id+2;
				}
				
				#pragma omp task
				ballAlg(l, l + (r - l) / 2, tree[id].left);
				#pragma omp task
				ballAlg(l + (r - l) / 2, r, tree[id].right);
			}
		}
	}
}

/* Prints the resulting tree to a file or stdout. */
void print_tree(node *tree) {
	FILE *fd = fopen("pts.txt", "w");
	fprintf(fd, "%d %ld\n", n_dims, 2 * np - 1);
	for (long i = 0; i < 2 * np - 1; ++i) {
		fprintf(fd, "%ld %ld %ld %f ", i, tree[i].left, tree[i].right, tree[i].radius);
		// fprintf(fd, "%f ", tree[i].radius);
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

	n_nodes = 2 * np - 1;
	// printf("n_nodes: %ld\n", n_nodes);
	n_levels = 1 + ceil(log(np)/log(2));
	// printf("n_levels: %ld\n", n_levels);
	index_last = 0;
	for (int i = 0; i < n_levels-1; i++)
		index_last += pow(2, i);
	// printf("index_last: %ld\n", index_last);
	n_last_level = n_nodes - index_last;
	// printf("n_last_level: %ld\n", n_last_level);
	l_lower = pow(2, n_levels-2) - 1;
	// printf("l_lower: %ld\n", l_lower);
	l_upper = index_last - 1;
	// printf("l_upper: %ld\n", l_upper);

	// omp_set_nested(1);

	#pragma omp parallel num_threads(4)
	{

		#pragma omp single
		{
			// Isto aqui me baixo não funciona como devia...
			int total_n_threads = omp_get_num_threads();
			printf("Total number of threads: %d\n", total_n_threads);

			// 1. Get input sample points
			pt_array = get_points(argc, argv, &n_dims, &np);

			n_nodes = 2 * np - 1;
			// printf("n_nodes: %ld\n", n_nodes);
			n_levels = 1 + ceil(log(np)/log(2));
			// printf("n_levels: %ld\n", n_levels);
			index_last = 0;
			for (int i = 0; i < n_levels-1; i++)
				index_last += pow(2, i);
			// printf("index_last: %ld\n", index_last);
			n_last_level = n_nodes - index_last;
			// printf("n_last_level: %ld\n", n_last_level);
			l_lower = pow(2, n_levels-2) - 1;
			// printf("l_lower: %ld\n", l_lower);
			l_upper = index_last - 1;
			// printf("l_upper: %ld\n", l_upper);

			#ifdef DEBUG
				FILE *fd = fopen("points.txt", "w");
				for (int i = 0; i < np; ++i) {
					print_point(pt_array[i], n_dims, stdout);
					print_point(pt_array[i], n_dims, fd);
				}
				fclose(fd);
			#endif

			idx = (long *)malloc(sizeof(long) * np);

			proj_scalar = (double *)malloc(np * sizeof(double));
			center1 = (double *)malloc(n_dims * sizeof(double));

			tree = (node *)malloc(n_nodes * sizeof(node));

		}
		
		#pragma omp for
			for (int i = 0; i < np; i++)
				idx[i] = i;

		#pragma omp for
			for (int i = 0; i < n_nodes; i++)
				tree[i].center = (double *)malloc(n_dims * sizeof(double));

	}

	ballAlg(0, np, 0);

	exec_time += omp_get_wtime();
	fprintf(stderr, "%.3lf\n", exec_time);

	print_tree(tree);

	free(center1);

	for (int i = 0; i < n_nodes; ++i)
		free(tree[i].center);
	free(tree);

	free(proj_scalar);
	free(idx);
	free(pt_array[0]);
	free(pt_array);

	return 0;
}
