#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>

#define RANGE 10
#define BLOCK_LOW(id,p,np) ((id)*(np)/(p))
#define BLOCK_HIGH(id,p,np) (BLOCK_LOW((id)+1,p,np)-1)
#define BLOCK_SIZE(id,p,np) (BLOCK_HIGH(id,p,np)-BLOCK_LOW(id,p,np)+1)

double ** pt_arr;
double * proj;
long * idx_local;
int id;
int p;
int p_initial;
long n_dims;
long np;
int * block_size;
double pivot;
int id_initial;
double m;

void print_point(double *pt, int dim) {
    for (int i = 0; i < dim; ++i) {
        printf("%f ", pt[i]);
    }
    printf("\n");
    fflush(stdout);
}

void swap(double *i, double *j) {
    double temp = *i;
    *i = *j;
    *j = temp;
}

int partition(double pivot, int l, int r, int piv_id) {

    printf("Before - id: %d ", id);
    print_point(proj, BLOCK_SIZE(id,p,np));

    if (r == l + 1) {
        if (proj[l] < pivot)
            return 1;
        return 0;
    }

    int i, j, low;
    if (id == piv_id) {
        i = l;
        low = l+1;
    }
    else {
        i = l-1;
        low = l;
    }
    for (j = low; j < r; j++)
    {
        if (proj[j] < pivot)
        {
            i++;
            swap(&proj[i], &proj[j]);
        }
    }

    if (id == piv_id)
        swap(&proj[i], &proj[0]);

    printf("After - i: %d, j: %d, id: %d ", i, j, id);
    print_point(proj, BLOCK_SIZE(id,p,np));

    return i+1;

}

void quick_select(int k, int l, int r, int p, int piv_id, MPI_Comm comm) {

    // MPI_Barrier(comm);
    // exit(0)

    printf("id: %d, l: %d, r: %d, k: %d, piv_id: %d\n", id, l, r, k, piv_id);
    fflush(stdout);

    int pivot = 0;
    while(piv_id < p) {
        if (id == piv_id) {
            if (r - l < 1)
                pivot = -1;
            else
                pivot = proj[l];
        }
        MPI_Bcast(&pivot, 1, MPI_DOUBLE, piv_id, comm);
        if (pivot == -1)
            piv_id++;
        else
            break;
    }

    sleep(1);

    int sum;
    int i = partition(pivot, l, r, piv_id);

    MPI_Allreduce(&i, &sum, 1, MPI_INT, MPI_SUM, comm);

    printf("sum: %d\n", sum);
    fflush(stdout);

    // sleep(1);

    if(sum == k){
        if (id == piv_id)
            m = proj[i-1];
        MPI_Bcast(&m, 1, MPI_DOUBLE, piv_id, comm);
        printf("id: %d, m: %f\n", id, m);
        // print_point(m, n_dims);
        fflush(stdout);      
    } else if(sum > k){
        quick_select(k, l, (id == piv_id) ? i-1 : i, p, piv_id, comm);
    } else {
        quick_select(k-sum, i, r, p, piv_id, comm);
    }

}

int main(int argc, char **argv) {

    MPI_Init (&argc, &argv);

    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);

    id_initial = id;
    p_initial = p;

    n_dims = atoi(argv[1]);
    np = atol(argv[2]);
    unsigned seed = atoi(argv[3]);
    srandom(seed);

    double * p_aux;
    // add extra space
    float extra = 5/4;
    p_aux = (double *)malloc(n_dims * int(ceil(np/p) * extra) * sizeof(double));
    pt_arr = (double **)malloc(int(ceil(np/p) * extra) * sizeof(double *));
    for (long i = 0; i < int(ceil(np/p) * extra); i++)
        pt_arr[i] = &p_aux[i * n_dims];

    block_size = (int *) malloc (p * sizeof(int));
    // branch_size = (int *) malloc (p * sizeof(int));
    idx_local = (long *) malloc (int(ceil(np/p) * extra)*sizeof(long));
    proj = (double *) malloc (int(ceil(np/p) * extra)*sizeof(double));
    // m = (double *) malloc (n_dims * sizeof(double));

    if (!id) {
        proj[0] = 45;
        proj[1] = 73;
        proj[2] = 16;
        proj[3] = 38;
    }
    else if (id == 1) {
        proj[0] = 58;
        proj[1] = 47;
        proj[2] = 22;
        proj[3] = 68;        
    }
    else if (id == 2) {
        proj[0] = 43;
        proj[1] = 19;
        proj[2] = 13;
        proj[3] = 24;        
    }
    else {
        proj[0] = 67;
        proj[1] = 34;
        proj[2] = 57;
        proj[3] = 23; 
    }

    int BL = BLOCK_LOW(id, p, np);
    int BH = BLOCK_HIGH(id, p, np);
    for (int i = 0; i < np; i++)
        for (int j = 0; j < n_dims; j++) {
            if (i >= BL && i <= BH) {
                pt_arr[i-BL][j] = RANGE * ((double)random()) / RAND_MAX;
                idx_local[i-BL] = i-BL;
            }
            else
                random();
        }

    // int previous = 0;
    for (int i = 0; i < p; i++) {
        block_size[i] = BLOCK_SIZE(i,p,np);
        // branch_size[i] = 2*block_size[i]-1 + previous;
        // previous = branch_size[i];
        // if(!id) {
        //     printf("i: %d, block_size[i]: %d, branch_size[i]: %d\n", i, block_size[i], branch_size[i]);
        //     fflush(stdout);
        // }
    }

    // for (int i = 0; i < BLOCK_SIZE(i,p,np); i++) {
    //     printf("id: %d ", id);
    //     fflush(stdout);
    //     print_point(pt_arr[i], n_dims);
    // }

    quick_select(7, 0, BLOCK_SIZE(id,p,np), p, 0, MPI_COMM_WORLD);

    free(block_size);
    // free(branch_size);
    free(idx_local);
    free(pt_arr[0]);
    free(pt_arr);

    MPI_Finalize();

}
