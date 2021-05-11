#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "gen_points.h"

#define BLOCK_LOW(id,p,np) ((id)*(np)/(p))
#define BLOCK_HIGH(id,p,np) (BLOCK_LOW((id)+1,p,np)-1)
#define BLOCK_SIZE(id,p,np) (BLOCK_HIGH(id,p,np)-BLOCK_LOW(id,p,np)+1)

int par_info = 0;

// This is used for the AllReduce
struct mpi_allreduce { 
    double max_d; 
    int idx; 
};

typedef struct tree_node {
    int center_idx;
    double radius;
    long left;
    long right;
    long node_id;
} node;

struct local_data {
    long local_idx;
    long global_idx;
    double proj;
};

long n_nodes; // Number of nodes
long n_levels; // Number of levels in the tree
long id_last; // First id of the last level
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
int tree_counter = 0;

MPI_Datatype mpi_data_struct; // MPI struct equivalent to local_data struct

// This is used for the scatterV of the sorted indexes
// (because we can't send it directly to the struct local_data)
long * idx_local;
// This is also used for the satterV
int * displs;
// This is a temp variable with the center for the processors to calculate the radius
// double * center;
// Another temp variable
int new_id;
// Aux variables for the all reduce
struct mpi_allreduce in, out;

void print_point(double *pt, int dim, FILE *fd) {
    for (int i = 0; i < dim; ++i)
        fprintf(fd, "%f ", pt[i]);
    fprintf(fd, "\n");
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

void print_struct(struct local_data * data, int size) {
    for (int i = 0; i < size; i++) {
        printf("id: %d, l_idx: %ld, g_idx: %ld, proj: %f\n", id, data[i].local_idx, data[i].global_idx, data[i].proj);
        fflush(stdout);
    }
}

void swap(long *i, long *j) {
    long temp = *i;
    *i = *j;
    *j = temp;
}

double dist(double *pt1, double *pt2) {
    double dist = 0.0;
    for (int d = 0; d < n_dims; ++d)
        dist += (pt1[d] - pt2[d]) * (pt1[d] - pt2[d]);
    return dist;
}

void furthest_apart(int l, int r, long *a_idx, long *b_idx) {
    double max_d = 0.0;
    double d;
    long idx0 = data[l].global_idx;
    for (int i = l + 1; i < r; ++i) {
        if (data[i].global_idx < idx0)
            idx0 = data[i].global_idx;
    }
    for (int i = l; i < r; ++i) {
        d = dist(pt_array[idx0], pt_array[data[i].global_idx]);
        if (d > max_d) {
            max_d = d;
            *a_idx = data[i].global_idx;
        }
    }
    max_d = 0;
    for (int i = l; i < r; ++i) {
        d = dist(pt_array[*a_idx], pt_array[data[i].global_idx]);
        if (d > max_d) {
            max_d = d;
            *b_idx = data[i].global_idx;
        }
    }
    if (pt_array[*a_idx][0] > pt_array[*b_idx][0])
        swap(a_idx, b_idx);
}

int get_kth_element(int k, int l, int h) {

    if (h == l + 1)
        return data[l].local_idx;

    else {
        int i = l;
        int j = h;

        while (i < j) {
            do {
                i++;
            } while (i < h && (data[i].proj <= data[l].proj));
            if (i == h) {
                j = h - 1;
                break;
            }
            do {
                j--;
            } while (data[j].proj > data[l].proj);
            if (i < j) {
                long g_idx_temp = data[i].global_idx;
                data[i].global_idx = data[j].global_idx;
                data[j].global_idx = g_idx_temp;
                double proj_temp = data[i].proj;
                data[i].proj = data[j].proj;
                data[j].proj = proj_temp;
            }
        }

        long g_idx_temp = data[l].global_idx;
        data[l].global_idx = data[j].global_idx;
        data[j].global_idx = g_idx_temp;
        double proj_temp = data[l].proj;
        data[l].proj = data[j].proj;
        data[j].proj = proj_temp;

        if (j - l == k)
            return data[j].local_idx;
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
        *median_2 = get_kth_element(0, l + (h - l) / 2, h);
    }
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

double orth_projv1(double *a, double *b, double *p) {
    double proj = 0.0;
    for (int i = 0; i < n_dims; ++i)
        proj += (p[i] - a[i]) * (b[i] - a[i]);

    return proj;
}

// Puts two sorted array in one sorted array
void merge_sort(struct local_data * data1, struct local_data * data2, struct local_data * data_merge, int size1, int size2) {

    int i = size1-1;
    int j = size2-1;
    int k = size1+size2;

    while(k > 0) {
        data_merge[--k] = (j < 0 || (i >= 0 && (data1[i].proj >= data2[j].proj))) ? data1[i--] : data2[j--];
    }

    // print_struct(data_merge, size1+size2);
}

// Puts sorted indexes into globalArray (defines only in process 0)
void mergeSort(int id, struct local_data * data, MPI_Comm comm, long * globalArray){
    
    int height = log(p) / log(2);
    int parent, rightChild, myHeight;
    struct local_data * data1;
    struct local_data * data2;
    struct local_data * data_merge;

    myHeight = 0;
    int size = get_block_size(id, myHeight);

    qsort(data, size, sizeof(struct local_data), compare);

    data1 = data;
    int new_size;

    while (myHeight < height) {
        parent = (id  & (~(1 << myHeight)));

        if (parent == id) {
            rightChild = (id | (1 << myHeight));

            data2 = (struct local_data *) malloc (size*sizeof(struct local_data));
            MPI_Recv(data2, get_block_size(rightChild, myHeight), mpi_data_struct, rightChild, 0, comm, MPI_STATUS_IGNORE);

            new_size = BLOCK_SIZE(id,p,np) + BLOCK_SIZE(rightChild,p,np);
            data_merge = (struct local_data *) malloc (new_size*sizeof(struct local_data));
            merge_sort(data1, data2, data_merge, size, get_block_size(rightChild, myHeight));
            data1 = data_merge;
            size = new_size;
            
            free(data2);
            data_merge = NULL;

            myHeight++;

        } else {

            MPI_Send(data1, get_block_size(id, myHeight), mpi_data_struct, parent, 0, comm);
            
            if(myHeight != 0) {
                free(data1);
            }
            myHeight = height;
        }
    }

    if(id == 0){
        for (int i = 0; i < np; i++) {
            globalArray[i] = data1[i].global_idx;
            if (par_info == 1) {
                printf("%ld ", globalArray[i]);
                fflush(stdout);                
            }

        }
        if (par_info == 1) {
            printf("\n");
            fflush(stdout);
        }
    }
}

void ballAlg_seq(long l, long r, long tree_id, int lvl) {

    // if (lvl == 2) {
    //     printf("SEQ - id: %d, l: %ld, r: %ld, tree_id: %ld, proc: %d\n", id, l, r, tree_id, p);
    //     fflush(stdout);
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     exit(0);
    // }

    // for (int i = l; i < r; ++i) {
    //     printf("Before - id: %d, l_idx: %ld, proj: %f\n", id, data[i].local_idx, data[i].proj);
    //     fflush(stdout);
    // }

    // MPI_Barrier(MPI_COMM_WORLD);
    // exit(0);

    if (r - l == 1) {
        tree[tree_counter].center_idx = l;
        tree[tree_counter].radius = 0;
        tree[tree_counter].left = -1;
        tree[tree_counter].right = -1;
        tree[tree_counter].node_id = tree_id;
        tree_counter++;
        return;
    }

    long a, b;
    // 2. Compute points a and b, furthest apart in the current set (approx)
    furthest_apart(l, r, &a, &b);

    // 3. Perform the orthogonal projection of all points onto line ab
    for (int i = l; i < r; ++i) {
        data[i].proj = orth_projv1(pt_array[a], pt_array[b], pt_array[data[i].global_idx]);
        printf("Before - id: %d, l_idx: %ld, proj: %f\n", id, data[i].local_idx, data[i].proj);
        fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // 4. Compute the center, defined as the median point over all projections
    int m1, m2 = -1;
    get_median(l, r, &m1, &m2);

    for (int i = l; i < r; ++i) {
        printf("After - id: %d, l_idx: %ld, proj: %f\n", id, data[i].local_idx, data[i].proj);
        fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    double abnorm = 0.0, aux, u;
    if ((r - l) % 2)
        u = orth_projv1(pt_array[a], pt_array[b], pt_array[data[m1].global_idx]);
    else
        u = (orth_projv1(pt_array[a], pt_array[b], pt_array[data[m1].global_idx]) + orth_projv1(pt_array[a], pt_array[b], pt_array[data[m2].global_idx])) / 2;
    for (int i = 0; i < n_dims; ++i) {
        aux = (pt_array[b][i] - pt_array[a][i]);
        centers[n_center][i] = u * aux;
        abnorm += aux * aux;
    }

    for (int i = 0; i < n_dims; ++i)
        centers[n_center][i] = centers[n_center][i] / abnorm + pt_array[a][i];

    double max_r = 0, rad;
    for (int i = l; i < r; ++i) {
        rad = dist(centers[n_center], pt_array[data[i].local_idx]);
        if (rad > max_r)
            max_r = rad;
    }

    if (lvl == 2) {
        printf("m1: %d, m2: %d\n", m1, m2);
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);
        exit(0);
    }

    else {
        tree[tree_counter].left = 2 * tree_id + 1;
        tree[tree_counter].right = 2 * tree_id + 2;
    }

    int temp_id = tree_counter;

    tree[tree_counter].node_id = tree_id;
    tree[tree_counter].center_idx = n_center;
    tree[tree_counter].radius = sqrt(max_r);

    fprintf(stdout, "%ld %ld %ld %f ", tree_id, tree[tree_counter].left, tree[tree_counter].right, tree[tree_counter].radius);
    print_point(centers[tree[tree_counter].center_idx], n_dims, stdout);

    n_center++;       
    tree_counter++;

    ballAlg_seq(l, l + (r - l) / 2, tree[temp_id].left, lvl + 1);
    ballAlg_seq(l + (r - l) / 2, r, tree[temp_id].right, lvl + 1);

}

double * center;
int p_initial;

void ballAlg(long l, long r, long tree_id, int lvl, int proc, MPI_Comm comm) {

    p = proc;

    // if (par_info == 1) {
    //     if (tree_id != 0) {
    //         printf("id: %d, l: %ld, r: %ld, tree_id: %ld, proc: %d\n", id, l, r, tree_id, p);
    //         fflush(stdout);
    //         MPI_Barrier(comm);
    //         // exit(0);
    //     }        
    // }

    if (!id) {

        if (r - l == 1) {
            tree[tree_counter].center_idx =  data[l].global_idx;
            tree[tree_counter].radius = 0;
            tree[tree_counter].left = -1;
            tree[tree_counter].right = -1;
            tree[tree_counter].node_id = tree_id;
            tree_counter++;
            return;
        }
    }

    // Local smaller index
    long idx0 = data[0].global_idx;
    // Global smaller index
    long idx0_global;

    for (int i = 0; i < BLOCK_SIZE(id,p,r-l); ++i) {
        if (data[i].global_idx < idx0)
            idx0 = data[i].global_idx;
    }
    
    MPI_Allreduce(&idx0, &idx0_global, 1, MPI_LONG, MPI_MIN, comm);

    long a, b;
    double d;
    // Local max distance
    in.max_d = 0;
    // Local index correspondent to the max distance
    in.idx = -1;
    for (int i = 0; i < BLOCK_SIZE(id,p,r-l); ++i) {
        d = dist(pt_array[idx0_global], pt_array[data[i].global_idx]);
        if (d > in.max_d) {
            in.max_d = d;
            in.idx = data[i].global_idx;
        }
    }

    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
    a = out.idx;

    in.max_d = 0;
    in.idx = -1;
    for (int i = 0; i < BLOCK_SIZE(id,p,r-l); ++i) {
        d = dist(pt_array[a], pt_array[data[i].global_idx]);
        if (d > in.max_d) {
            in.max_d = d;
            in.idx = data[i].global_idx;
        }
    }

    // MPI_Bcast(center, n_dims, MPI_DOUBLE, 0, comm);
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
    
    b = out.idx;

    if (pt_array[a][0] > pt_array[b][0])
        swap(&a, &b);

    if (tree_id != 0) {
        if (par_info == 1) {
            printf("node_id: %ld, id: %d, a: %ld, b: %ld\n", tree_id, id, a, b);
            fflush(stdout);
        }
    }

    for (int i = 0; i < BLOCK_SIZE(id,p,r-l); ++i) {
        data[i].proj = orth_projv1(pt_array[a], pt_array[b], pt_array[data[i].global_idx]);
        if (tree_id != 0) {
            if (par_info == 1) {
                printf("node_id: %ld, id: %d, idx[i]: %ld, proj[i]: %f, i: %d\n", tree_id, id, data[i].global_idx, data[i].proj, i);
                fflush(stdout);
            }
        }
    }

    if (!id) {
        mergeSort(id, data, comm, idx_sorted);
    }
    else {
        mergeSort(id, data, comm, NULL);
    }

    if (!id) {
        if (par_info == 1) {
            printf("node_id: %ld - ", tree_id);
            for (int i = 0; i < np; i++) {
                printf("%ld ", idx_sorted[i]);
                fflush(stdout);
            }
            printf("\n");
            fflush(stdout);
        }

        int m1, m2 = -1;
        if ((r - l) % 2 == 1) 
            m1 = idx_sorted[(r - l - 1) / 2];
        else {
            m1 = idx_sorted[(r - l) / 2 - 1];
            m2 = idx_sorted[(r - l) / 2];
        }

        if (par_info == 1) {
            printf("n_center: %ld, m1: %d, m2: %d\n", n_center, m1, m2);
            fflush(stdout);
        }

        double abnorm = 0.0, aux, u;
        if ((r - l) % 2)
            u = orth_projv1(pt_array[a], pt_array[b], pt_array[m1]);
        else
            u = (orth_projv1(pt_array[a], pt_array[b], pt_array[m1]) + orth_projv1(pt_array[a], pt_array[b], pt_array[m2])) / 2;
        for (int i = 0; i < n_dims; ++i) {
            aux = (pt_array[b][i] - pt_array[a][i]);
            centers[n_center][i] = u * aux;
            abnorm += aux * aux;
        }

        for (int i = 0; i < n_dims; ++i)
            centers[n_center][i] = centers[n_center][i] / abnorm + pt_array[a][i];
        center = centers[n_center];
        if (par_info == 1) {
            printf("center: %f %f\n", center[0], center[1]);
            fflush(stdout);
        }
    }

    // // ---
    // MPI_Barrier(MPI_COMM_WORLD);
    // double v1 = 1;
    // double v2 = 0;
    // MPI_Reduce(&v1, &v2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // printf("id: %d, v2: %f\n", id, v2);
    // MPI_Barrier(MPI_COMM_WORLD);
    // //

    MPI_Bcast(center, n_dims, MPI_DOUBLE, 0, comm);

    // printf("%ld\n", BLOCK_SIZE(id,p,r-l));
    // fflush(stdout);
    // MPI_Barrier(comm);
    // printf("%f %f\n", center[0], center[1]);
    // fflush(stdout);
    // MPI_Barrier(comm);

    if (par_info == 1) {
        for (int i = 0; i < BLOCK_SIZE(id,p,r-l); i++) {
            printf("id: %d, dist: %f\n", id, dist(center, pt_array[data[i].global_idx]));
            fflush(stdout);
        }
        MPI_Barrier(comm);
    }

    double max_d = 0;
    for (int i = 0; i < BLOCK_SIZE(id,p,r-l); i++) {
        d = dist(center, pt_array[data[i].global_idx]);
        if (d > max_d)
            max_d = d;
    }

    double rad = 0;

    if (par_info == 1) {
        printf("id: %d, max_d: %f, rad: %f\n", id, max_d, rad);
        fflush(stdout);
        MPI_Barrier(comm);
    }
   
    MPI_Reduce(&max_d, &rad, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

    if (par_info == 1) {
        printf("id: %d, rad: %f\n", id, rad);
        fflush(stdout);
    }

    if (((np & (np - 1)) != 0) && lvl == n_levels - 1) {
        tree[tree_counter].left = id_last;
        tree[tree_counter].right = id_last + 1;
        id_last += 2;
    } else {
        tree[tree_counter].left = 2 * tree_id + 1;
        tree[tree_counter].right = 2 * tree_id + 2;
    }

    int temp_id = tree_counter;

    if (!id) {

        tree[tree_counter].node_id = tree_id;
        tree[tree_counter].center_idx = n_center;
        tree[tree_counter].radius = sqrt(rad);
        n_center++;

        // if (par_info == 1) {
            fprintf(stdout, "%ld %ld %ld %f ", tree_id, tree[tree_counter].left, tree[tree_counter].right, tree[tree_counter].radius);
            print_point(centers[tree[tree_counter].center_idx], n_dims, stdout);
        // }
                
        tree_counter++;
    }

    MPI_Scatterv(idx_sorted, block_size, displs, MPI_LONG, idx_local, block_size[id], MPI_LONG, 0, comm);

    for (int i = 0; i < BLOCK_SIZE(id,p,r-l); i++) {
        data[i].local_idx = i;
        data[i].global_idx = idx_local[i];
        printf("id: %d, tree_id: %ld, l_idx: %ld, g_idx: %ld\n", id, tree_id, data[i].local_idx, data[i].global_idx);
        fflush(stdout);
    }

    // // HELP
    // // Fill the block size
    // for (int i = 0; i < new_proc; i++) {
    //     block_size[i] = BLOCK_SIZE(i,new_proc,r-l);
    //     displs[i] = l + BLOCK_LOW(i,new_proc,r-l);
    // }

    if (lvl == log(p_initial)/log(2)-1) {
        if (!id)
            ballAlg_seq(0, BLOCK_SIZE(id,p,r-l), tree[temp_id].left, lvl + 1);
        else
            ballAlg_seq(0, BLOCK_SIZE(id,p,r-l), tree[temp_id].right, lvl + 1);
    } else {
        MPI_Comm new_comm;

        MPI_Comm_split(comm, id/2, id, &new_comm);

        if (tree_id != 0)
            MPI_Comm_free(&comm);

        MPI_Comm_rank(new_comm, &new_id);
        MPI_Comm_size(new_comm, &p);

        if (id/2 < 1) {
            id = new_id;
            printf("id: %d, l: %ld, r: %ld, tree_id: %ld, proc: %d\n", id, l, r, tree_id, p);
            fflush(stdout);
            ballAlg(l, l + (r - l) / 2, tree[temp_id].left, lvl + 1, p, new_comm);
        }
        else {
            id = new_id;
            printf("id: %d, l: %ld, r: %ld, tree_id: %ld, proc: %d\n", id, l, r, tree_id, p);
            fflush(stdout);
            ballAlg(l + (r - l) / 2, r, tree[temp_id].right, lvl + 1, p, new_comm);
        }
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

    p_initial = p;

    // // ---

    // MPI_Barrier(MPI_COMM_WORLD);

    // double v1 = 1;
    // double v2 = 0;

    // MPI_Reduce(&v1, &v2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // printf("id: %d, v2: %f\n", id, v2);

    // MPI_Barrier(MPI_COMM_WORLD);

    // // --- 

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
    // Array of structs with indices and projections
    data = (struct local_data *) malloc (ceil(np/p)*sizeof(struct local_data));
    idx_local = (long *) malloc (ceil(np/p)*sizeof(long));
    displs = (int *) malloc (ceil(np/p)*sizeof(int));

    // Fill the block size
    for (int i = 0; i < p; i++) {
        block_size[i] = BLOCK_SIZE(i,p,np);
        displs[i] = BLOCK_LOW(i,p,np);
    }

    // Fill the struct array
    for (int i = 0; i < BLOCK_SIZE(id,p,np); ++i) {
        data[i].local_idx = i;
        data[i].global_idx = BLOCK_LOW(id,p,np) + i;
    }

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

    // See if this will be needed or not    
    n_nodes = 2 * np - 1;
    n_levels = ceil(log(np) / log(2));
    id_last = pow(2, n_levels) - 1;

    // Alocation of part of the tree in all processors
    tree = (node *)malloc((log(p)/log(2)+(2*(np/pow(2, log(p)/log(2)-1))-1)) * sizeof(node));
    double *center_aux = (double *)malloc((np - 1) * n_dims * sizeof(double));

    // Alocation of the center in all processors
    centers = (double **)malloc((log(p)/log(2)+np/pow(2, log(p)/log(2))-1-(np/p)) * sizeof(double *));
    for (int i = 0; i < ceil(np/p); ++i)
        centers[i] = &center_aux[i * n_dims];
    // Alocation of the sorted indexes in all processors
    idx_sorted = (long *) malloc (np * sizeof(long));
    // Temporary center for radius calculation
    center = (double *) malloc (n_dims * sizeof(double));

    ballAlg(0, np, 0, 0, p, MPI_COMM_WORLD);

    elapsed_time += MPI_Wtime();

    MPI_Type_free(&mpi_data_struct);
    MPI_Finalize ();

    fprintf(stderr, "%.1lf\n", elapsed_time);

    // print_tree(tree);

    free(block_size);
    free(data);
    free(idx_local);
    free(displs);
    free(idx_sorted);
    free(center);

    free(pt_array[0]);
    free(pt_array);

    free(centers[0]);
    free(centers);

    free(tree);

    exit(0);
}
