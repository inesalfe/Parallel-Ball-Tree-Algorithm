#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "gen_points.h"

#define BLOCK_LOW(id,p,np) ((id)*(np)/(p))
#define BLOCK_HIGH(id,p,np) (BLOCK_LOW((id)+1,p,np)-1)
#define BLOCK_SIZE(id,p,np) (BLOCK_HIGH(id,p,np)-BLOCK_LOW(id,p,np)+1)

struct { 
    double max_d; 
    int idx; 
} in, out;

typedef struct tree_node {
    int center_idx;
    double radius;
    long left;
    long right;
} node;

struct local_data {
    long local_idx;
    long global_idx;
    double proj;
};

long n_nodes; // Number of nodes
long n_levels; // Number of levels in the tree
long id_last; // First id of the last level
long c_id; // Id of the center
int id, p; // Id of the processor and number of processors
double **centers; // Array of centers
long n_center = 0; // 
int n_dims; // Number of dimensions
long np; // Number of points
double **pt_array; // Array of points
node *tree; // Array of tree nodes
struct local_data * data; // Array of structs local to each processor
long * idx_sorted; // Global array of sorted global indexes
int * block_size; // Array of the number of elements in each processor

MPI_Datatype mpi_data_struct; // MPI struct equivalent to local_data struct

void print_struct(struct local_data * data, int size) {
    for (int i = 0; i < size; i++) {
        printf("l_idx: %ld, g_idx: %ld, proj: %f\n", data[i].local_idx, data[i].global_idx, data[i].proj);
        fflush(stdout);
    }
}

void swap(long *i, long *j) {
    long temp = *i;
    *i = *j;
    *j = temp;
}

int compare (const void * a, const void * b) {
   return ( (*(struct local_data*)a).proj - (*(struct local_data*)b).proj );
}

int get_block_size(int id, int cur_height) {
    int size = 0;
    for (int i = id; i < id + pow(2, cur_height); i++) {
        size += block_size[i];
    }
    return size;
}

void print_point(double *pt, int dim, FILE *fd) {
    for (int i = 0; i < dim; ++i)
        fprintf(fd, "%f ", pt[i]);
    fprintf(fd, "\n");
}

double dist(double *pt1, double *pt2) {
    double dist = 0.0;
    for (int d = 0; d < n_dims; ++d)
        dist += (pt1[d] - pt2[d]) * (pt1[d] - pt2[d]);
    return dist;
}

double orth_projv1(double *a, double *b, double *p) {
    double proj = 0.0;
    for (int i = 0; i < n_dims; ++i)
        proj += (p[i] - a[i]) * (b[i] - a[i]);

    return proj;
}

void merge_sort(struct local_data * data1, struct local_data * data2, struct local_data * data_merge, int size1, int size2) {

    int i = size1-1;
    int j = size2-1;
    int k = size1+size2;

    while(k > 0) {
        data_merge[--k] = (j < 0 || (i >= 0 && (data1[i].proj >= data2[j].proj))) ? data1[i--] : data2[j--];
    }

    print_struct(data_merge, size1+size2);
}

void mergeSort(int id, struct local_data * data, MPI_Comm comm, long * globalArray){
    
    int height = log(p) / log(2);
    int parent, rightChild, myHeight;
    struct local_data * data1;
    struct local_data * data2;
    struct local_data * data_merge;

    myHeight = 0;
    int size = get_block_size(id, myHeight);

    qsort(data, size, sizeof(struct local_data), compare);
    print_struct(data, size);

    data1 = data;
    int new_size;

    while (myHeight < height) {
        parent = (id & (~(1 << myHeight)));

        if (parent == id) {
            rightChild = (id | (1 << myHeight));

            data2 = (struct local_data *) malloc (size*sizeof(struct local_data));
            MPI_Recv(data2, get_block_size(rightChild, myHeight), mpi_data_struct, rightChild, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            new_size = BLOCK_SIZE(id,p,np) + BLOCK_SIZE(rightChild,p,np);
            data_merge = (struct local_data *) malloc (new_size*sizeof(struct local_data));
            merge_sort(data1, data2, data_merge, size, get_block_size(rightChild, myHeight));
            data1 = data_merge;
            size = new_size;
            
            free(data2);
            data_merge = NULL;

            myHeight++;

        } else {

            MPI_Send(data1, get_block_size(id, myHeight), mpi_data_struct, parent, 0, MPI_COMM_WORLD);
            
            if(myHeight != 0) {
                free(data1);
            }
            myHeight = height;
        }
    }

    if(id == 0){
        for (int i = 0; i < np; i++) {
            globalArray[i] = data1[i].global_idx;
            printf("%ld ", globalArray[i]);
            fflush(stdout);
        }
        printf("\n");
        fflush(stdout);
    }
}

void ballAlg(long l, long r, long tree_id, int lvl) {

    if (!id) {

        if (r - l == 1) {
            tree[tree_id].center_idx = l;
            tree[tree_id].radius = 0;
            tree[tree_id].left = -1;
            tree[tree_id].right = -1;
            return;
        }

        c_id = n_center++;
    }

    long idx0 = data[0].global_idx;
    long idx0_global;

    for (int i = 0; i < BLOCK_SIZE(id,p,np); ++i) {
        if (data[i].global_idx < idx0)
            idx0 = data[i].global_idx;
    }
    
    MPI_Allreduce(&idx0, &idx0_global, 1, MPI_LONG, MPI_MIN, MPI_COMM_WORLD);

    long a, b;
    double d;
    in.max_d = 0;
    in.idx = -1;
    for (int i = 0; i < BLOCK_SIZE(id,p,np); ++i) {
        d = dist(pt_array[idx0_global], pt_array[data[i].global_idx]);
        if (d > in.max_d) {
            in.max_d = d;
            in.idx = data[i].global_idx;
        }
    }

    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    a = out.idx;

    in.max_d = 0;
    in.idx = -1;
    for (int i = 0; i < BLOCK_SIZE(id,p,np); ++i) {
        d = dist(pt_array[a], pt_array[data[i].global_idx]);
        if (d > in.max_d) {
            in.max_d = d;
            in.idx = data[i].global_idx;
        }
    }

    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

    b = out.idx;

    if (pt_array[a][0] > pt_array[b][0])
        swap(&a, &b);
    
    printf("id: %d, a: %ld, b: %ld\n", id, a, b);
    fflush(stdout);

    for (int i = 0; i < BLOCK_SIZE(id,p,np); ++i) {
        data[i].proj = orth_projv1(pt_array[a], pt_array[b], pt_array[data[i].global_idx]);
        printf("id: %d, idx[i]: %ld, proj[i]: %f, i: %d\n", id, data[i].global_idx, data[i].proj, i);
        fflush(stdout);
    }

    if (id == 0) {
        idx_sorted = (long *) malloc (np * sizeof(long));
        mergeSort(id, data, MPI_COMM_WORLD, idx_sorted);
    }
    else {
        mergeSort(id, data, MPI_COMM_WORLD, NULL);
    }

    if (id == 0) {
        for (int i = 0; i < np; i++) {
            printf("%ld ", idx_sorted[i]);
            fflush(stdout);
        }
        printf("\n");
        fflush(stdout);
        free(idx_sorted);
    }

    double * center = (double *) malloc (n_dims * sizeof(double));

    if (!id) {
        int m1, m2 = -1;
        if ((r - l) % 2 == 1) 
            m1 = idx_sorted[(r - l - 1) / 2];
        else {
            m1 = idx_sorted[(r - l) / 2 - 1];
            m2 = idx_sorted[(r - l) / 2];
        }        
        double abnorm = 0.0, aux, u;
        if ((r - l) % 2)
            u = orth_projv1(pt_array[a], pt_array[b], pt_array[m1]);
        else
            u = (orth_projv1(pt_array[a], pt_array[b], pt_array[m1]) + orth_projv1(pt_array[a], pt_array[b], pt_array[m2])) / 2;
        for (int i = 0; i < n_dims; ++i) {
            aux = (pt_array[b][i] - pt_array[a][i]);
            centers[c_id][i] = u * aux;
            abnorm += aux * aux;
        }
        for (int i = 0; i < n_dims; ++i) {
            centers[c_id][i] = centers[c_id][i] / abnorm + pt_array[a][i];
            printf("%f ", centers[c_id][i]);
            fflush(stdout);
        }
        printf("\n");
        center = centers[c_id];
    }

    MPI_Bcast(center, n_dims, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double max_d = 0;
    for (int i = 0; i < BLOCK_SIZE(id,p,np); ++i) {
        d = dist(center, pt_array[data[i].global_idx]);
        if (d > max_d)
            max_d = d;
    }
    
    double rad;
    MPI_Reduce(&max_d, &rad, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (!id) {

        tree[tree_id].center_idx = c_id;
        tree[tree_id].radius = sqrt(rad);

        if (((np & (np - 1)) != 0) && lvl == n_levels - 1) {
            tree[tree_id].left = id_last;
            tree[tree_id].right = id_last + 1;
            id_last += 2;
        } else {
            tree[tree_id].left = 2 * tree_id + 1;
            tree[tree_id].right = 2 * tree_id + 2;
        }

        fprintf(stdout, "%ld %ld %ld %f ", tree_id, tree[tree_id].left, tree[tree_id].right, tree[tree_id].radius);
        print_point(centers[tree[tree_id].center_idx], n_dims, stdout);

        MPI_Barrier(MPI_COMM_WORLD);
        exit(0);

        ballAlg(l, l + (r - l) / 2, tree[tree_id].left, lvl + 1);
        ballAlg(l + (r - l) / 2, r, tree[tree_id].right, lvl + 1);
    
    }

}

void print_tree(node *tree) {
    fprintf(stdout, "%d %ld\n", n_dims, 2 * np - 1);
    for (long i = 0; i < 2 * np - 1; ++i) {
        fprintf(stdout, "%ld %ld %ld %f ", i, tree[i].left, tree[i].right, tree[i].radius);
        if (tree[i].left == -1)
            print_point(pt_array[tree[i].center_idx], n_dims, stdout);
        else
            print_point(centers[tree[i].center_idx], n_dims, stdout);
    }
}

int main(int argc, char **argv) {

    // Value for time counting
    double elapsed_time;

    // Initialization od MPI
    MPI_Init (&argc, &argv);

    // Time starts counting here
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    // Get id of processor
    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    // Get number of processors
    MPI_Comm_size (MPI_COMM_WORLD, &p);

    // Check if the arguments are corrects
    if (argc != 4) {
        if (!id)
            printf ("Command line: %s <m>\n", argv[0]);
        MPI_Finalize();
        exit (1);
    }

    // All processors have acess to all the points
    pt_array = get_points(argc, argv, &n_dims, &np);

    // Every processor knows how much data the other processors have
    // This is needed for the merge sort
    block_size = (int *) malloc (p * sizeof(int));
    for (int i = 0; i < p; i++)
        block_size[i] = BLOCK_SIZE(i,p,np);

    // Creating of a MPI Struct to pass data around
    const int nitems = 3;
    int blocklengths[3] = {1,1,1};
    MPI_Datatype types[3] = {MPI_LONG, MPI_LONG, MPI_DOUBLE};
    MPI_Aint offsets[3];

    offsets[0] = offsetof(struct local_data, local_idx);
    offsets[1] = offsetof(struct local_data, global_idx);
    offsets[2] = offsetof(struct local_data, proj);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_data_struct);
    MPI_Type_commit(&mpi_data_struct);

    // Array of structs with indices and projections
    data = (struct local_data *) malloc (ceil(np/p)*sizeof(struct local_data));

    // Fill the struct array
    for (int i = 0; i < BLOCK_SIZE(id,p,np); ++i) {
        data[i].local_idx = i;
        data[i].global_idx = BLOCK_LOW(id,p,np)+i;
        printf("id: %d, pt_array[0]: %f, pt_array[1]: %f, idx: %ld\n", id, pt_array[data[i].global_idx][0], pt_array[data[i].global_idx][1], data[i].global_idx);
        fflush(stdout);
    }

    // See if this will be needed or not    
    n_nodes = 2 * np - 1;
    n_levels = ceil(log(np) / log(2));
    id_last = pow(2, n_levels) - 1;

    // Alocation of the tree in all processors
    tree = (node *)malloc((2 * np - 1) * sizeof(node));
    double *center_aux = (double *)malloc((np - 1) * n_dims * sizeof(double));

    // The center will the calculated only in the first processor
    if (!id) {
        centers = (double **)malloc((np - 1) * sizeof(double *));
        for (int i = 0; i < np - 1; ++i)
            centers[i] = &center_aux[i * n_dims];
    }

    ballAlg(0, np, 0, 0);

    elapsed_time += MPI_Wtime();

    MPI_Type_free(&mpi_data_struct);
    MPI_Finalize ();

    fprintf(stderr, "%.1lf\n", elapsed_time);

    print_tree(tree);

    free(pt_array[0]);
    free(pt_array);

    free(centers[0]);
    free(centers);

    free(tree);

    exit(0);
}
