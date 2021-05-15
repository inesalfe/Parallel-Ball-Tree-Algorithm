#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>

struct reduce_d_i { 
    double max_d; 
    int idx; 
} in_max, out_max;

struct reduce_i_i { 
    int idx0; 
    int p_id; 
} in_min, out_min;

long n_dims;
long np;
int p;
int p_initial;
int id;
int id_initial;
int aux;
double * p_aux;
double ** pt_arr;
int * block_size;
long * idx_global;
double * pt_arr_idx0;
double * pt_arr_a;
double * pt_arr_b;
double * proj;
double * pivots;
double * glob_pivots;
int * send_displs;
int * recv_displs;
int * send_counts;
int * recv_counts;
double * recv_buffer;
double * p_aux_2;
double * center;
double **centers;
long n_center = 0;

#define RANGE 10
#define BLOCK_LOW(id,p,np) ((id)*(np)/(p))
#define BLOCK_HIGH(id,p,np) (BLOCK_LOW((id)+1,p,np)-1)
#define BLOCK_SIZE(id,p,np) (BLOCK_HIGH(id,p,np)-BLOCK_LOW(id,p,np)+1)

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

double orth_projv1(double *a, double *b, double *p) {
    double proj = 0.0;
    for (int i = 0; i < n_dims; ++i)
        proj += (p[i] - a[i]) * (b[i] - a[i]);
    return proj;
}

void ballAlg(long tree_id, int lvl, MPI_Comm comm) {

    if (!id) {

        if (block_size[id_initial] == 1) {
            // Não sei o que fazer no center_idx
            // tree[tree_counter].center_idx = ;
            // tree[tree_counter].radius = 0;
            // tree[tree_counter].left = -1;
            // tree[tree_counter].right = -1;
            // tree[tree_counter].node_id = tree_id;
            // tree_counter++;
            return;
        }
    }

    in_min.idx0 = idx_global[0];
    for (int i = 0; i < block_size[id_initial]; i++) {
        if (idx_global[i] < in_min.idx0) {
            in_min.idx0 = idx_global[i];
            in_min.p_id = id;
        }
    }

    MPI_Allreduce(&in_min, &out_min, 1, MPI_2INT, MPI_MINLOC, comm);

    int root_id;
    if (id == out_min.p_id) {
        root_id = id;
        // printf("idx0: %d ", out_min.p_id);
        // fflush(stdout);
        // print_point(pt_arr[out_min.idx0], n_dims);
        for (int i = 0; i < n_dims; i++) {
            pt_arr_idx0[i] = pt_arr[out_min.idx0][i];
        }
    }

    MPI_Bcast(pt_arr_idx0, n_dims, MPI_DOUBLE, root_id, comm);

    // printf("id: %d ", id);
    // fflush(stdout);
    // print_point(pt_arr_idx0, n_dims);

    double d;
    in_max.max_d = 0;
    in_max.idx = -1;
    for (int i = 0; i < block_size[id_initial]; ++i) {
        d = dist(pt_arr_idx0, pt_arr[i]);
        if (d > in_max.max_d) {
            in_max.max_d = d;
            in_max.idx = i;
        }
    }

    // Rever a partir daqui - é precio os outros processadores saberem qual é o root

    MPI_Allreduce(&in_max, &out_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
    
    if (in_max.max_d == out_max.max_d) {
        root_id = id;
        for (int i = 0; i < n_dims; i++) {
            pt_arr_a[i] = pt_arr[out_max.idx][i];
        }
        printf("id: %d ", id);
        fflush(stdout);
        print_point(pt_arr_a, n_dims);
    }

    MPI_Bcast(pt_arr_a, n_dims, MPI_DOUBLE, root_id, comm);

    printf("id: %d ", id);
    fflush(stdout);
    print_point(pt_arr_a, n_dims);

    in_max.max_d = 0;
    in_max.idx = -1;
    for (int i = 0; i < block_size[id_initial]; ++i) {
        d = dist(pt_arr_a, pt_arr[i]);
        if (d > in_max.max_d) {
            in_max.max_d = d;
            in_max.idx = i;
        }
    }

    MPI_Allreduce(&in_max, &out_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
    
    if (in_max.max_d == out_max.max_d)    
        for (int i = 0; i < n_dims; i++) {
            pt_arr_b[i] = pt_arr[out_max.idx][i];
        }

    MPI_Bcast(pt_arr_b, n_dims, MPI_DOUBLE, 0, comm);

    if (pt_arr_a[0] > pt_arr_b[0]) {
        double temp;
        for (int i = 0; i < n_dims; i++) {
            temp = pt_arr_a[i];
            pt_arr_a[i] = pt_arr_a[i];
            pt_arr_b[i] = pt_arr_a[i];
        }
    }

    for (int i = 0; i < block_size[id_initial]; ++i) {
        proj[i] = orth_projv1(pt_arr_a, pt_arr_b, pt_arr[i]);
    }

    int min_index = 0;
    for(int i = 0; i < block_size[id_initial]; ++i) {
        if(proj[i] < proj[min_index]){
            min_index = i;
        }
    }
    swap(0, min_index);
    pivots[0] = proj[0];

    int temp_idx = -1;
    int piv_idx, prev_piv_idx = 0;
    for(int i = 1; i < p; ++i) {
        piv_idx = i * block_size[id_initial] / p - prev_piv_idx - 1;
        temp_idx = quickselect_seq(piv_idx, prev_piv_idx+1, block_size[id_initial]);
        prev_piv_idx = temp_idx;
        pivots[i] = proj[temp_idx];
    }

    MPI_Gather(pivots, p, MPI_DOUBLE, glob_pivots, p, MPI_DOUBLE, 0, comm);

    if(!id){
        prev_piv_idx = -1;
        for(int i = 1; i < p; ++i) {
            piv_idx = i * p - prev_piv_idx - 1;
            pivots[i-1] = quickselect_seq_2(piv_idx, prev_piv_idx+1, p*p);
            prev_piv_idx += piv_idx + 1;
        }
    }

    MPI_Bcast(pivots, p-1, MPI_DOUBLE, 0, comm);

    send_counts[0] = partition(pivots[0], 0, block_size[id_initial]);
    int sum = send_counts[0];
    for(int i = 1; i < p-1; i++) {
        send_counts[i] = partition(pivots[i], sum, block_size[id_initial]) - sum;
        sum += send_counts[i];
    }
    send_counts[p-1] = block_size[id_initial] - sum;

    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, comm);

    int sum1 = 0;
    int sum2 = 0;
    for(int i = 0; i < p; i++) {
        send_displs[i] = sum2;
        recv_displs[i] = sum1;
        sum1 += recv_counts[i];
        sum2 += send_counts[i];
    }

    MPI_Alltoallv(proj, send_counts, send_displs, MPI_DOUBLE, recv_buffer, recv_counts, recv_displs, MPI_DOUBLE, comm);

    free(proj);
    proj = recv_buffer;
    
    for(int i = 0; i < p; ++i){
        send_counts[i] *= n_dims;
        send_displs[i] *= n_dims;
        recv_counts[i] *= n_dims;
        recv_displs[i] *= n_dims;
    }

    MPI_Alltoallv(p_aux, send_counts, send_displs, MPI_DOUBLE, p_aux_2, recv_counts, recv_displs, MPI_DOUBLE, comm);

    for (long i = 0; i < sum1; i++)
        pt_arr[i] = &p_aux_2[i * n_dims];

    int size = (id < p/2 + (p%2)) ? sum1 : 0;
    int total_size = 0;

    printf("id: %d, size: %d, total_size: %d\n", id, size, total_size);
    fflush(stdout);

    MPI_Allreduce(&size, &total_size, 1, MPI_INT, MPI_SUM, comm);

    MPI_Barrier(MPI_COMM_WORLD);
    exit(0);

    int m1 = -1, m2 = -1;
    int min_proj;
    double u;
    if ((id == p/2) && (p%2 == 1)) {
        m1 = quickselect_seq(np/2-!(np%2)-(total_size-sum1), 0, sum1);
        // printf("id: %d, (np)/2-(total_size-sum1): %ld\n", id, np/2-!(np%2)-(total_size-sum1));
        if (np%2 == 0) {
            min_proj = proj[m1+1];
            m2 = m1+1;
            for (int i = m1+1; block_size[id_initial]; i++)
                if (proj[i] < min_proj) {
                    min_proj = proj[i];
                    m2 = i;
                }
            swap(m2, m1+1);
            u = (proj[m1] + proj[m1+1])/2;
            MPI_Send(proj+m1+1, sum1-(m1+1), MPI_DOUBLE, id+1, 0, comm);
            MPI_Send(p_aux_2+(m1+1)*n_dims, (sum1-(m1+1))*n_dims, MPI_DOUBLE, id+1, 1, comm);
            sum1 = m1+1;
        }
        else {
            u = proj[m1];
            // printf("id: %d, u: %f\n", id, u);
            // fflush(stoutd);
            MPI_Send(proj+m1, sum1-m1, MPI_DOUBLE, id+1, 0, comm);
            MPI_Send(p_aux_2+m1*n_dims, (sum1-m1)*n_dims, MPI_DOUBLE, id+1, 1, comm);
            sum1 = m1;
        }
        // printf("id: %d, m1: %d, m2: %d, proj[m1]: %f\n", id, m1, m2, proj[m1]);
    }
    else if ((id == p/2+1) && (p%2 == 1)) {
        u = 0;
        MPI_Recv(proj+sum1, total_size-(np/2), MPI_DOUBLE, id-1, 0, comm, MPI_STATUS_IGNORE);        
        MPI_Recv(p_aux_2+sum1*n_dims, (total_size-(np/2))*n_dims, MPI_DOUBLE, id-1, 1, comm, MPI_STATUS_IGNORE);
        sum1 += total_size-(np/2);
        for (long i = 0; i < sum1; i++) {
            pt_arr[i] = &p_aux_2[i * n_dims];
        }
    }
    else if ((np%2 == 1) && (p%2 == 0)) {
        if (id == p/2 - 1) {
            if (total_size >= np/2) {
                m1 = quickselect_seq((np)/2-(total_size-sum1), 0, sum1);
                u = proj[m1];
                MPI_Send(proj+m1, sum1-m1, MPI_DOUBLE, id+1, 0, comm);
                MPI_Send(p_aux_2+m1*n_dims, (sum1-m1)*n_dims, MPI_DOUBLE, id+1, 1, comm);
                sum1 = m1;
            }
            else {
                u = 0;
                fflush(stdout);
                MPI_Recv(proj+sum1, (np/2)-total_size, MPI_DOUBLE, id+1, 0, comm, MPI_STATUS_IGNORE);        
                MPI_Recv(p_aux_2+sum1*n_dims, ((np/2)-total_size)*n_dims, MPI_DOUBLE, id-1, 1, comm, MPI_STATUS_IGNORE);
                sum1 += (np/2)-total_size;
                for (long i = 0; i < sum1; i++) {
                    pt_arr[i] = &p_aux_2[i * n_dims];
                }
            }
        }
        else if (id == p/2) {
            if (total_size < np/2) {
                m1 = quickselect_seq((np)/2-(total_size), 0, sum1);
                u = proj[m1];
                MPI_Send(proj, m1, MPI_DOUBLE, id-1, 0, comm);
                MPI_Send(p_aux_2, m1*n_dims, MPI_DOUBLE, id-1, 1, comm);
                sum1 = sum1-m1;
                for (int i = 0; i < sum1; i++) {
                    proj[i] = proj[i+m1+1];
                    // pt_arr[i] = &p_aux_2[(i+m1+1) * n_dims];
                    pt_arr[i] = pt_arr[i+m1+1];
                }
            }
            else {
                u = 0;
                MPI_Recv(proj+sum1, total_size-(np/2), MPI_DOUBLE, id+1, 0, comm, MPI_STATUS_IGNORE);        
                MPI_Recv(p_aux_2+sum1*n_dims, (total_size-(np/2))*n_dims, MPI_DOUBLE, id+1, 1, comm, MPI_STATUS_IGNORE);
                sum1 += total_size-(np/2);
                for (long i = 0; i < sum1; i++) {
                    pt_arr[i] = &p_aux_2[i * n_dims];
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
                    for (int i = m1+1; i < block_size[id_initial]; i++)
                        if (proj[i] < min_proj) {
                            min_proj = proj[i];
                            m2 = i;
                        }
                    swap(m2, m1+1);
                    u = (proj[m1]+proj[m1+1])/2;
                    MPI_Send(proj+m2, sum1-m2, MPI_DOUBLE, id+1, 0, comm);
                    MPI_Send(p_aux_2+m2*n_dims, (sum1-m2)*n_dims, MPI_DOUBLE, id+1, 1, comm);
                    sum1 = sum1-m2;
                }
                else {
                    u = proj[m1]/2;
                }
            }
            else {
                u = 0;
                MPI_Recv(proj+sum1, (np/2)-total_size, MPI_DOUBLE, id+1, 0, comm, MPI_STATUS_IGNORE);        
                MPI_Recv(p_aux_2+sum1*n_dims, ((np/2)-total_size)*n_dims, MPI_DOUBLE, id+1, 1, comm, MPI_STATUS_IGNORE);
                sum1 += (np/2)-total_size;
                for (long i = 0; i < sum1; i++) {
                    pt_arr[i] = &p_aux_2[i * n_dims];
                }
            }
        }
        else if (id == p/2) {
            if (total_size < np/2) {
                m1 = quickselect_seq(np/2-1-total_size, 0, sum1);
                min_proj = proj[m1+1];
                m2 = m1+1;
                for (int i = m1+1; i < block_size[id_initial]; i++)
                    if (proj[i] < min_proj) {
                        min_proj = proj[i];
                        m2 = i;
                    }
                swap(m2, m1+1);
                u = (proj[m1]+proj[m1+1])/2;
                MPI_Send(proj, m1+1, MPI_DOUBLE, id-1, 0, comm);
                MPI_Send(p_aux_2, (m1+1)*n_dims, MPI_DOUBLE, id-1, 1, comm);
                sum1 = sum1-(m1+1);
                for (int i = 0; i < sum1; i++) {
                    proj[i] = proj[i+m1+1];
                    // pt_arr[i] = &p_aux_2[(i+m1+1) * n_dims];
                    pt_arr[i] = pt_arr[i+m1+1];
                }
            }
            else if (total_size == np/2) {
                m1 = -1;
                min_proj = proj[m1+1];
                m2 = m1+1;
                for (int i = m1+1; i < block_size[id_initial]; i++)
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
                MPI_Recv(p_aux_2+sum1*n_dims, (total_size-(np/2))*n_dims, MPI_DOUBLE, id-1, 1, comm, MPI_STATUS_IGNORE);
                sum1 += total_size-(np/2);
                for (long i = 0; i < sum1; i++) {
                    pt_arr[i] = &p_aux_2[i * n_dims];
                }
            }
        }
    }
    else {
        u = 0;
    }

    block_size[id_initial] = sum1;

    double u_aux;
    MPI_Allreduce(&u, &u_aux, 1, MPI_DOUBLE, MPI_SUM, comm);
    u = u_aux;

    double abnorm = 0.0;
    for (int i = 0; i < n_dims; ++i) {
        u_aux = (pt_arr_b[i] - pt_arr_a[i]);
        center[i] = u * u_aux;
        abnorm += u_aux * u_aux;
    }

    for (int i = 0; i < n_dims; ++i) {
        center[i] = center[i] / abnorm + pt_arr_a[i];
    }

    double max_d = 0;
    for (int i = 0; i < block_size[id_initial]; i++) {
        d = dist(center, pt_arr[i]);
        if (d > max_d)
            max_d = d;
    }

    double rad = 0;
   
    MPI_Reduce(&max_d, &rad, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

    MPI_Barrier(MPI_COMM_WORLD);
    exit(0);
    
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

    n_dims = atoi(argv[1]);
    np = atol(argv[2]);
    unsigned seed = atoi(argv[3]);
    srandom(seed);

    p_initial = p;
    id_initial = id;

    aux = ceil(np/p) * 5/2;
    p_aux = (double *)malloc(n_dims * aux * sizeof(double));
    pt_arr = (double **)malloc(aux * sizeof(double *));
    for (long i = 0; i < aux; i++)
        pt_arr[i] = &p_aux[i * n_dims];

    block_size = (int *) malloc (p * sizeof(int));
    for (int i = 0; i < p; i++) {
        block_size[i] = BLOCK_SIZE(i,p,np);
        // printf("id: %d, i: %d, block_size[i]: %d\n", id, i, block_size[i]);
        // fflush(stdout);
    }

    idx_global = (long *) malloc (aux * sizeof(long));

    int BL = BLOCK_LOW(id, p, np);
    int BH = BLOCK_HIGH(id, p, np);
    for (int i = 0; i < np; i++)
        for (int j = 0; j < n_dims; j++) {
            if (i >= BL && i <= BH) {
                pt_arr[i-BL][j] = RANGE * ((double)random()) / RAND_MAX;
            }
            else
                random();
        }

    // Se desse para englobar este no ciclo de cima era fixe
    for (int i = 0; i < block_size[id]; i++) {
        idx_global[i] = BLOCK_LOW(id,p,np) + i;
        // printf("id: %d, i: %d, idx_global[i]: %ld\n", id, i, idx_global[i]);
        // fflush(stdout);
    }

    pt_arr_idx0 = (double *) malloc (n_dims*sizeof(double));
    pt_arr_a = (double *) malloc (n_dims*sizeof(double));
    pt_arr_b = (double *) malloc (n_dims*sizeof(double));
    proj = (double *) malloc (aux * sizeof(double));
    pivots = (double *) malloc (p * sizeof(double));
    glob_pivots = (double *) malloc (p * p * sizeof(double));
    send_displs = (int *) malloc (p * sizeof(int));
    recv_displs = (int *) malloc (p * sizeof(int));
    send_counts = (int *) malloc (p * sizeof(int));
    recv_counts = (int *) malloc (p * sizeof(int));
    recv_buffer = (double *) malloc (aux * sizeof(double));
    p_aux_2 = (double *)malloc(n_dims * aux * sizeof(double));

    // Fazer contas para ver quantos centros alocar (eu meti só o tamanho máximo)
    double *center_aux = (double *)malloc((np - 1) * n_dims * sizeof(double));
    centers = (double **)malloc((np-1) * sizeof(double *));
    for (int i = 0; i < np-1; ++i)
        centers[i] = &center_aux[i * n_dims];
    center = (double *) malloc (n_dims * sizeof(double));

    ballAlg(0, 0, MPI_COMM_WORLD);

    elapsed_time += MPI_Wtime();

    free(p_aux);
    free(pt_arr);
    free(block_size);
    free(idx_global);
    free(pt_arr_idx0);
    free(pt_arr_a);
    free(pt_arr_b);
    free(proj);
    free(pivots);
    free(glob_pivots);
    free(send_displs);
    free(recv_displs);
    free(send_counts);
    free(recv_counts);
    free(recv_buffer);
    free(p_aux_2);

    if (!id_initial) {
        fprintf(stderr, "%.1lf\n", elapsed_time);
        //print_tree(tree);
    }

    MPI_Finalize();

}