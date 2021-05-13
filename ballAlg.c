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
double ** pt_arr_2;
double * proj_2;
double * p_aux_2;
long * idx_local;
int id;
int id_initial;
int p;
int p_initial;
long n_dims;
long np;
int * block_size;
int aux;
double * pivots;
double * glob_pivots;
int * send_displs;
int * recv_displs;

void print_point(double *pt, int dim) {
    for (int i = 0; i < dim; ++i) {
        printf("%f ", pt[i]);
    }
    printf("\n");
    fflush(stdout);
}

void swap(int i, int j) {
    double temp = proj[i];
    proj[i] = proj[j];
    proj[j] = temp;
    for (int k = 0; k < n_dims; k++) {
        double temp_pt = pt_arr[i][k];
        pt_arr[i][k] = pt_arr[j][k];
        pt_arr[j][k] = temp_pt;
    }
}

void swap_vec(int i, int j) {
    double temp = glob_pivots[i];
    glob_pivots[i] = glob_pivots[j];
    glob_pivots[j] = temp;
}

int partition(double pivot, int l, int r) {

    if (r == l + 1) {
        if (proj[l] <= pivot)
            return 1;
        else
            return 0;
    }

    int i = l-1;
    int j;
    for (j = l; j < r; j++)
    {
        if (proj[j] < pivot)
        {
            i++;
            swap(i, j);
        }
    }

    if (proj[i+1] == pivot)
        return i+2;
    else
        return i+1;

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

double quickselect_seq(int k, int l, int h) {

    if (h == l + 1)
        return proj[l];

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
            return proj[j];
        if (k < j - l)
            return quickselect_seq(k, l, j);
        else
            return quickselect_seq(k - (j - l) - 1, j + 1, h);
    }
}

void get_median(int l, int r, int p, MPI_Comm comm){
    
    printf("id: %d, l: %d, r: %d, p: %d\n", id, l, r, p);
    fflush(stdout);

    // Ecnontrar primeiro pivot - distância mínima
    int min_index = l;
    int i;
    for(i = l; i < r; ++i) {
        if(proj[i] < proj[min_index]){
            min_index = i;
        }
    }
    swap(0, min_index);
    pivots[0] = proj[0];

    printf("id: %d, pivots[0]: %f\n", id, pivots[0]);
    fflush(stdout);

    // MPI_Barrier(MPI_COMM_WORLD);
    // exit(0);

    // Outros dois pivots são escolhido como o primeiro elemento dado um índice
    for(int i = 1; i < p; ++i) {
        pivots[i] = quickselect_seq((r-l)/p-1, (i-1)*(r-l)/p+1, r);
    }

    // print_point(proj, 9);

    // printf("id: %d, pivots[0]: %f, pivots[1]: %f, pivots[2]: %f\n", id, pivots[0], pivots[1], pivots[2]);
    // fflush(stdout);

    MPI_Barrier(MPI_COMM_WORLD);
    // exit(0);

    // Enviar os pivots locais para o primeiro processador
    MPI_Gather(pivots, p, MPI_DOUBLE, glob_pivots, p, MPI_DOUBLE, 0, comm);

    // if (!id) {
    //     print_point(glob_pivots, p*p);
    // }

    // O primeiro processador calcula p-1 pivots utilizando quick select
    if(!id){
        for(int i = 1; i < p; ++i) {
            pivots[i-1] = quickselect_seq_2((r-l)/p-1, (i-1)*(r-l)/p+1, r);
        }
        // print_point(glob_pivots, 9);
    }

    // printf("id: %d, pivots[0]: %f, pivots[1]: %f\n", id, pivots[0], pivots[1]);
    // fflush(stdout);

    // Enviam-se esses pivots para todos os processadores
    MPI_Bcast(pivots, p-1, MPI_DOUBLE, 0, comm);

    // print_point(proj, 9);  

    // Partições segundos os pivots
    pivots[0] = partition(pivots[0], 0, r);
    int sum = pivots[0];
    for(int i = 1; i < p-1; i++) {
        pivots[i] = partition(pivots[i], (int)pivots[i-1]+1, r) - pivots[i-1];
        sum += pivots[i];
    }
    pivots[p-1] = r - sum;

    printf("id: %d, pivots[0]: %f, pivots[1]: %f, pivots[2]: %f\n", id, pivots[0], pivots[1], pivots[2]);
    fflush(stdout);
    // print_point(proj, 9);    

    // Array of sizes
    MPI_Alltoall(pivots, 1, MPI_DOUBLE, glob_pivots, 1, MPI_DOUBLE, comm);

    print_point(glob_pivots, 3);

    int sum1 = 0;
    int sum2 = 0;
    for(int i = 0; i < p; i++) {
        send_displs[i] = sum2;
        recv_displs[i] = sum1;
        sum1 += glob_pivots[i];
        sum2 += pivots[i];
    }

    printf("id: %d, sum: %d\n", id, sum);

    // Fazer free disto to fim desta função
    if (sum > aux) {
        p_aux_2 = (double *)malloc(n_dims * (sum-aux) * sizeof(double));
        pt_arr_2 = (double **)malloc((sum-aux) * sizeof(double *));
        for (long i = 0; i < sum-aux; i++)
            pt_arr_2[i] = &p_aux_2[i * n_dims];
        proj_2 = (double *) malloc ((sum-aux) * sizeof(double));
    }
    else {
        MPI_Alltoallv(proj, pivots, send_displs, MPI_DOUBLE, proj, glob_pivots, recv_displs, MPI_DOUBLE, comm);
   
    }

    MPI_Barrier(MPI_COMM_WORLD);
    exit(0);

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

    float extra = 5/4;
    aux = ceil(np/p) * extra;
    double * p_aux = (double *)malloc(n_dims * aux * sizeof(double));
    pt_arr = (double **)malloc(aux * sizeof(double *));
    for (long i = 0; i < aux; i++)
        pt_arr[i] = &p_aux[i * n_dims];

    block_size = (int *) malloc (p * sizeof(int));
    pivots = (double *) malloc (p * sizeof(double));
    glob_pivots = (double *) malloc (p * p * sizeof(double));
    idx_local = (long *) malloc (aux * sizeof(long));
    proj = (double *) malloc (aux * sizeof(double));
    send_displs = (int *) malloc (p * sizeof(int));
    recv_displs = (int *) malloc (p * sizeof(int));

    if (!id) {
        proj[0] = 15;
        proj[1] = 46;
        proj[2] = 48;
        proj[3] = 6;
        proj[4] = 93;
        proj[5] = 39;
        proj[6] = 72;
        proj[7] = 91;
        proj[8] = 14;
    }
    else if (id == 1) {
        proj[0] = 36;
        proj[1] = 69;
        proj[2] = 40;
        proj[3] = 89;
        proj[4] = 61;
        proj[5] = 97;
        proj[6] = 12;
        proj[7] = 21;
        proj[8] = 54;       
    }
    else {
        proj[0] = 53;
        proj[1] = 97;
        proj[2] = 84;
        proj[3] = 58;
        proj[4] = 32;
        proj[5] = 27;
        proj[6] = 33;
        proj[7] = 72;
        proj[8] = 20;         
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

    // for (int i = 0; i < BLOCK_SIZE(i,p,np); i++)
    //     print_point(pt_arr[i], n_dims);

    for (int i = 0; i < p; i++) {
        block_size[i] = BLOCK_SIZE(i,p,np);
    }

    get_median(0, BLOCK_SIZE(id,p,np), p, MPI_COMM_WORLD);

    free(block_size);
    free(idx_local);
    free(pivots);
    free(glob_pivots);
    free(proj);
    free(send_displs);
    free(recv_displs);
    free(pt_arr[0]);
    free(pt_arr);

    MPI_Finalize();

}
