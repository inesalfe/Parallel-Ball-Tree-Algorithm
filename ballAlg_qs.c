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
double * recv_buffer;
double * p_aux_3;
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
int * send_counts;
int * recv_counts;
double * p_aux;

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
    double temp_pt;
    for (int k = 0; k < n_dims; k++) {
        temp_pt = pt_arr[i][k];
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
        if (proj[j] <= pivot)
        {
            i++;
            swap(i, j);
        }
    }

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

    // printf("id: %d, pivots[0]: %f\n", id, pivots[0]);
    // fflush(stdout);

    MPI_Barrier(MPI_COMM_WORLD);
    // exit(0);

    // Outros dois pivots são escolhido como o primeiro elemento dado um índice
    int idx = -1;
    int piv_idx, prev_piv_idx = 0;
    for(int i = 1; i < p; ++i) {
        piv_idx = i * (r-l) / p - prev_piv_idx - 1;
        idx = quickselect_seq(piv_idx, prev_piv_idx+1, r);
        prev_piv_idx = idx;
        pivots[i] = proj[idx];
    }

    // if (!id)
    //     print_point(proj, 6);  
    // else
    //     print_point(proj, 7);

    // printf("id: %d, pivots[0]: %f, pivots[1]: %f, pivots[2]: %f, pivots[3]: %f\n", id, pivots[0], pivots[1], pivots[2], pivots[3]);
    // fflush(stdout);

    MPI_Barrier(MPI_COMM_WORLD);

    // Enviar os pivots locais para o primeiro processador
    MPI_Gather(pivots, p, MPI_DOUBLE, glob_pivots, p, MPI_DOUBLE, 0, comm);

    // if (!id) {
    //     print_point(glob_pivots, p*p);
    // }

    // O primeiro processador calcula p-1 pivots utilizando quick select
    if(!id){
        prev_piv_idx = -1;
        for(int i = 1; i < p; ++i) {
            piv_idx = i * p - prev_piv_idx - 1;
            pivots[i-1] = quickselect_seq_2(piv_idx, prev_piv_idx+1, p*p);
            // printf("pivots[i-1]: %f, piv_idx: %d, prev_piv_idx: %d\n", pivots[i-1], piv_idx, prev_piv_idx);
            prev_piv_idx += piv_idx + 1;
        }
        // print_point(glob_pivots, 16);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    // exit(0);

    // printf("id: %d, pivots[0]: %f, pivots[1]: %f\n", id, pivots[0], pivots[1]);
    // fflush(stdout);

    // Enviam-se esses pivots para todos os processadores
    MPI_Bcast(pivots, p-1, MPI_DOUBLE, 0, comm);

    // if (!id)
    //     print_point(proj, 6);  
    // else
    //     print_point(proj, 7);

    // Partições segundos os pivots
    send_counts[0] = partition(pivots[0], 0, r);
    int sum = send_counts[0];
    for(int i = 1; i < p-1; i++) {
        send_counts[i] = partition(pivots[i], sum, r) - sum;
        sum += send_counts[i];
    }

    send_counts[p-1] = r - sum;

    // printf("id: %d, send_counts[0]: %d, send_counts[1]: %d, send_counts[2]: %d, send_counts[3]: %d\n", id, send_counts[0], send_counts[1], send_counts[2], send_counts[3]);
    // fflush(stdout);

    // if (!id)
    //     print_point(proj, 6);  
    // else
    //     print_point(proj, 7);

    // print_point(proj, 9);    

    // Array of sizes
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, comm);

    // printf("id: %d, recv_counts[0]: %d, recv_counts[1]: %d, recv_counts[2]: %d, recv_counts[3]: %d\n", id, recv_counts[0], recv_counts[1], recv_counts[2], recv_counts[3]);
    // fflush(stdout);

    int sum1 = 0;
    int sum2 = 0;
    for(int i = 0; i < p; i++) {
        send_displs[i] = sum2;
        recv_displs[i] = sum1;
        sum1 += recv_counts[i];
        sum2 += send_counts[i];
    }

    printf("id: %d, sum1: %d\n", id, sum1);

    MPI_Barrier(MPI_COMM_WORLD);
    // exit(0);

    // printf("id: %d, recv_displs[0]: %d, recv_displs[1]: %d, recv_displs[2]: %d, recv_displs[3]: %d\n", id, recv_displs[0], recv_displs[1], recv_displs[2], recv_displs[3]);
    // fflush(stdout);
    // printf("id: %d, send_displs[0]: %d, send_displs[1]: %d, send_displs[2]: %d, send_displs[3]: %d\n", id, send_displs[0], send_displs[1], send_displs[2], send_displs[3]);
    // fflush(stdout);
    // printf("aux: %d\n", aux);
    // fflush(stdout);

    //Fazer free disto to fim desta função
    // if (sum1 > aux) {
    //     p_aux_2 = (double *)malloc(n_dims * (sum1-aux) * sizeof(double));
    //     pt_arr_2 = (double **)malloc((sum1-aux) * sizeof(double *));
    //     for (long i = 0; i < sum1-aux; i++)
    //         pt_arr_2[i] = &p_aux_2[i * n_dims];
    //     proj_2 = (double *) malloc ((sum1-aux) * sizeof(double));
    // }

    MPI_Alltoallv(proj, send_counts, send_displs, MPI_DOUBLE, recv_buffer, recv_counts, recv_displs, MPI_DOUBLE, comm);
    // print_point(recv_buffer, sum1); 

    free(proj);
    proj = recv_buffer;

    for(int i = 0; i < p; ++i){
        send_counts[i] *= n_dims;
        send_displs[i] *= n_dims;
        recv_counts[i] *= n_dims;
        recv_displs[i] *= n_dims;
    }

    // printf("id: %d, recv_displs[0]: %d, recv_displs[1]: %d, recv_displs[2]: %d\n", id, recv_displs[0], recv_displs[1], recv_displs[2]);
    // printf("id: %d, send_displs[0]: %d, send_displs[1]: %d, send_displs[2]: %d\n", id, send_displs[0], send_displs[1], send_displs[2]);
    // printf("id: %d, recv_counts[0]: %d, recv_counts[1]: %d, recv_counts[2]: %d\n", id, recv_counts[0], recv_counts[1], recv_counts[2]);
    // printf("id: %d, send_counts[0]: %d, send_counts[1]: %d, send_counts[2]: %d\n", id, send_counts[0], send_counts[1], send_counts[2]);

    // for (int i = 0; i < BLOCK_SIZE(i,p,np); i++)
    //     print_point(pt_arr[i], n_dims);

    MPI_Alltoallv(p_aux, send_counts, send_displs, MPI_DOUBLE, p_aux_3, recv_counts, recv_displs, MPI_DOUBLE, comm);
    for (long i = 0; i < sum1; i++)
        pt_arr[i] = &p_aux_3[i * n_dims];
    // for (int i = 0; i < sum1; i++)
    //     print_point(p_aux_3+i*n_dims, n_dims);

    int size = (id < p/2 + (p%2)) ? sum1 : 0;
    int total_size = 0;
    MPI_Allreduce(&size, &total_size, 1, MPI_INT, MPI_SUM, comm);

    printf("t_size: %d\n", total_size);

    // MPI_Barrier(MPI_COMM_WORLD);
    // exit(0);

    // print_point(proj, sum1);

    int m1 = -1, m2 = -1;
    int min_proj;
    double u;
    if ((id == p/2) && (p%2 == 1)) {
        m1 = quickselect_seq(np/2-!(np%2)-(total_size-sum1), 0, sum1);
        // printf("id: %d, (np)/2-(total_size-sum1): %ld\n", id, np/2-!(np%2)-(total_size-sum1));
        if (np%2 == 0) {
            min_proj = proj[m1+1];
            m2 = m1+1;
            for (int i = m1+1; i < r; i++)
                if (proj[i] < min_proj) {
                    min_proj = proj[i];
                    m2 = i;
                }
            swap(m2, m1+1);
            u = (proj[m1] + proj[m1+1])/2;
            // Enviar de m1+1 a r para a direita
            MPI_Send(proj+m1+1, sum1-(m1+1), MPI_DOUBLE, id+1, 0, comm);
            MPI_Send(p_aux_3+(m1+1)*n_dims, (sum1-(m1+1))*n_dims, MPI_DOUBLE, id+1, 1, comm);
            sum1 = m1+1;
        }
        else {
            u = proj[m1];
            // printf("id: %d, u: %f\n", id, u);
            // fflush(stdout);
            // Enviar de m1 a r para a direita
            MPI_Send(proj+m1, sum1-m1, MPI_DOUBLE, id+1, 0, comm);
            MPI_Send(p_aux_3+m1*n_dims, (sum1-m1)*n_dims, MPI_DOUBLE, id+1, 1, comm);
            sum1 = m1;
        }
        // printf("id: %d, m1: %d, m2: %d, proj[m1]: %f\n", id, m1, m2, proj[m1]);
    }
    else if ((id == p/2+1) && (p%2 == 1)) {
        u = 0;
        MPI_Recv(proj+sum1, total_size-(np/2), MPI_DOUBLE, id-1, 0, comm, MPI_STATUS_IGNORE);        
        MPI_Recv(p_aux_3+sum1*n_dims, (total_size-(np/2))*n_dims, MPI_DOUBLE, id-1, 1, comm, MPI_STATUS_IGNORE);
        sum1 += total_size-(np/2);
        for (long i = 0; i < sum1; i++) {
            pt_arr[i] = &p_aux_3[i * n_dims];
        }
    }
    else if ((np%2 == 1) && (p%2 == 0)) {
        if (id == p/2 - 1) {
            if (total_size >= np/2) {
                m1 = quickselect_seq((np)/2-(total_size-sum1), 0, sum1);
                u = proj[m1];
                MPI_Send(proj+m1, sum1-m1, MPI_DOUBLE, id+1, 0, comm);
                MPI_Send(p_aux_3+m1*n_dims, (sum1-m1)*n_dims, MPI_DOUBLE, id+1, 1, comm);
                sum1 = m1;
            }
            else {
                u = 0;
                fflush(stdout);
                MPI_Recv(proj+sum1, (np/2)-total_size, MPI_DOUBLE, id+1, 0, comm, MPI_STATUS_IGNORE);        
                MPI_Recv(p_aux_3+sum1*n_dims, ((np/2)-total_size)*n_dims, MPI_DOUBLE, id-1, 1, comm, MPI_STATUS_IGNORE);
                sum1 += (np/2)-total_size;
                for (long i = 0; i < sum1; i++) {
                    pt_arr[i] = &p_aux_3[i * n_dims];
                }
            }
        }
        else if (id == p/2) {
            if (total_size < np/2) {
                m1 = quickselect_seq((np)/2-(total_size), 0, sum1);
                u = proj[m1];
                MPI_Send(proj, m1, MPI_DOUBLE, id-1, 0, comm);
                MPI_Send(p_aux_3, m1*n_dims, MPI_DOUBLE, id-1, 1, comm);
                sum1 = sum1-m1;
                for (int i = 0; i < sum1; i++) {
                    proj[i] = proj[i+m1+1];
                    // pt_arr[i] = &p_aux_3[(i+m1+1) * n_dims];
                    pt_arr[i] = pt_arr[i+m1+1];
                }
            }
            else {
                u = 0;
                MPI_Recv(proj+sum1, total_size-(np/2), MPI_DOUBLE, id+1, 0, comm, MPI_STATUS_IGNORE);        
                MPI_Recv(p_aux_3+sum1*n_dims, (total_size-(np/2))*n_dims, MPI_DOUBLE, id+1, 1, comm, MPI_STATUS_IGNORE);
                sum1 += total_size-(np/2);
                for (long i = 0; i < sum1; i++) {
                    pt_arr[i] = &p_aux_3[i * n_dims];
                }
            }
        }
        else {
            u = 0;
        }
    }
    else if ((np%2 == 0) && (p%2 == 0)) {
        if (id == p/2 - 1) {
            if (total_size >= np/2) {
                m1 = quickselect_seq((np)/2-1-(total_size-sum1), 0, sum1);
                if (total_size != np/2) {
                    min_proj = proj[m1+1];
                    m2 = m1+1;
                    for (int i = m1+1; i < r; i++)
                        if (proj[i] < min_proj) {
                            min_proj = proj[i];
                            m2 = i;
                        }
                    swap(m2, m1+1);
                    u = (proj[m1]+proj[m1+1])/2;
                    MPI_Send(proj+m2, sum1-m2, MPI_DOUBLE, id+1, 0, comm);
                    MPI_Send(p_aux_3+m2*n_dims, (sum1-m2)*n_dims, MPI_DOUBLE, id+1, 1, comm);
                    sum1 = sum1-m2;
                }
                else {
                    u = proj[m1]/2;
                }
            }
            else {
                u = 0;
                MPI_Recv(proj+sum1, (np/2)-total_size, MPI_DOUBLE, id+1, 0, comm, MPI_STATUS_IGNORE);        
                MPI_Recv(p_aux_3+sum1*n_dims, ((np/2)-total_size)*n_dims, MPI_DOUBLE, id+1, 1, comm, MPI_STATUS_IGNORE);
                sum1 += (np/2)-total_size;
                for (long i = 0; i < sum1; i++) {
                    pt_arr[i] = &p_aux_3[i * n_dims];
                }
            }
        }
        else if (id == p/2) {
            if (total_size < np/2) {
                m1 = quickselect_seq(np/2-1-total_size, 0, sum1);
                min_proj = proj[m1+1];
                m2 = m1+1;
                for (int i = m1+1; i < r; i++)
                    if (proj[i] < min_proj) {
                        min_proj = proj[i];
                        m2 = i;
                    }
                swap(m2, m1+1);
                u = (proj[m1]+proj[m1+1])/2;
                MPI_Send(proj, m1+1, MPI_DOUBLE, id-1, 0, comm);
                MPI_Send(p_aux_3, (m1+1)*n_dims, MPI_DOUBLE, id-1, 1, comm);
                sum1 = sum1-(m1+1);
                for (int i = 0; i < sum1; i++) {
                    proj[i] = proj[i+m1+1];
                    // pt_arr[i] = &p_aux_3[(i+m1+1) * n_dims];
                    pt_arr[i] = pt_arr[i+m1+1];
                }
            }
            else if (total_size == np/2) {
                m1 = -1;
                min_proj = proj[m1+1];
                m2 = m1+1;
                for (int i = m1+1; i < r; i++)
                    if (proj[i] < min_proj) {
                        min_proj = proj[i];
                        m2 = i;
                    }
                swap(m2, m1+1);
                u = proj[0]/2;
            }
            else {
                u = 0;
                MPI_Recv(proj+sum1, total_size-(np/2), MPI_DOUBLE, id-1, 0, comm, MPI_STATUS_IGNORE);        
                MPI_Recv(p_aux_3+sum1*n_dims, (total_size-(np/2))*n_dims, MPI_DOUBLE, id-1, 1, comm, MPI_STATUS_IGNORE);
                sum1 += total_size-(np/2);
                for (long i = 0; i < sum1; i++) {
                    pt_arr[i] = &p_aux_3[i * n_dims];
                }
            }
        }
    }
    else {
        u = 0;
    }

    print_point(proj, sum1);

    for (int i = 0; i < sum1; i++)
        print_point(pt_arr[i], n_dims);

    if(id == 1) {
        printf("id: %d, u: %f\n", id, u);
        fflush(stdout);
    }

    double u_aux;
    MPI_Allreduce(&u, &u_aux, 1, MPI_DOUBLE, MPI_SUM, comm);
    u = u_aux;

    printf("id: %d, m1: %d, m2: %d, u: %f\n", id, m1, m2, u);

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

    aux = ceil(np/p) * 5/2;
    p_aux = (double *)malloc(n_dims * aux * sizeof(double));
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
    send_counts = (int *) malloc (p * sizeof(int));
    recv_counts = (int *) malloc (p * sizeof(int));
    recv_buffer = (double *) malloc (aux * sizeof(double));
    p_aux_3 = (double *)malloc(n_dims * aux * sizeof(double));

    if (!id) {
        proj[0] = 15;
        proj[1] = 46;
        proj[2] = 48;
        proj[3] = 6;
        proj[4] = 93;
        proj[5] = 39;
    }
    else if (id == 1) {
        proj[0] = 72;
        proj[1] = 91;
        proj[2] = 14;
        proj[3] = 36;
        proj[4] = 69;
        proj[5] = 40;
        proj[6] = 89;
    }
    else if (id == 2) {
        proj[0] = 61;
        proj[1] = 97;
        proj[2] = 12;
        proj[3] = 21;
        proj[4] = 54;           
        proj[5] = 53;
        proj[6] = 97;
    }
    else if (id == 3) {
        proj[0] = 84;
        proj[1] = 58;
        proj[2] = 32;
        proj[3] = 27;
        proj[4] = 33;
        proj[5] = 72;
        proj[6] = 20;
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
    free(send_counts);
    free(recv_counts);
    free(send_displs);
    free(recv_displs);
    free(pt_arr[0]);
    free(pt_arr);

    MPI_Finalize();

}
