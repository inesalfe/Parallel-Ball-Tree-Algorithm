#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
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

long n_dims;
long np;
int p;
int p_initial;
int id;
int id_initial;
int new_id;
int aux;
double * p_aux;
int * block_size;
int * branch_size;
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
long * recv_buffer_long;
double * p_aux_2;
double * center;
double ** centers;
long n_center = 0;
long tree_counter = 0;
node * tree;
pt *pts;

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
        temp_pt = p_aux[i*n_dims+k];
        p_aux[i*n_dims+k] = p_aux[j*n_dims+k];
        p_aux[j*n_dims+k] = temp_pt;
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
            return l+1;
        else
            return l;
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

long tree_idx;

long ballAlg_seq(long l, long r) {

    // printf("id_initial: %d, tree_idx: %ld, l: %ld, r: %ld, tree_counter: %ld, n_center: %ld\n", id_initial, tree_idx, l, r, tree_counter, n_center);
    // fflush(stdout);

    // MPI_Barrier(MPI_COMM_WORLD);
    // exit(0);

    if (r - l == 1) {
        tree[tree_counter].center_idx = l;
        tree[tree_counter].radius = 0;
        tree[tree_counter].left = -1;
        tree[tree_counter].right = -1;
        tree[tree_counter].node_id = tree_idx;
        // printf("%ld %ld %ld %f ", tree[tree_counter].node_id, tree[tree_counter].left, tree[tree_counter].right, tree[tree_counter].radius);
        // print_point(pts[l].pt, n_dims);
        // fflush(stdout);
        printf("%f ", tree[tree_counter].radius);
        print_point(pts[l].pt, n_dims);
        fflush(stdout);
        tree_idx++;
        tree_counter++;
        return tree_idx-1;
    }

    long temp_counter = tree_counter++;

    long a, b;
    furthest_apart(l, r, &a, &b);

    // MPI_Barrier(MPI_COMM_WORLD);
    // exit(0);

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

    tree[temp_counter].center_idx = n_center;
    tree[temp_counter].radius = sqrt(max_r);
    tree[temp_counter].node_id = tree_idx;

    // printf("%ld %f ", tree[temp_counter].node_id, tree[temp_counter].radius);
    // fflush(stdout);
    // print_point(centers[tree[temp_counter].center_idx], n_dims);

    printf("%f ", tree[temp_counter].radius);
    fflush(stdout);
    print_point(centers[tree[temp_counter].center_idx], n_dims);

    tree_idx++;
    n_center++;

    tree[temp_counter].left = ballAlg_seq(l, l + (r - l) / 2);
    tree[temp_counter].right = ballAlg_seq(l + (r - l) / 2, r);
    // printf("%ld %ld %ld %f ", tree[tree_counter].node_id, tree[tree_counter].left, tree[tree_counter].right, tree[tree_counter].radius);
    // print_point(centers[tree[tree_counter].center_idx], n_dims);
    // fflush(stdout);
}

void ballAlg(long n_points, long tree_id, int lvl, MPI_Comm comm) {

    // if (!id) {
    //     if (block_size[id_initial] == 1) {
    //         tree[tree_counter].center_idx = ;
    //         tree[tree_counter].radius = 0;
    //         tree[tree_counter].left = -1;
    //         tree[tree_counter].right = -1;
    //         tree[tree_counter].node_id = tree_id;
    //         tree_counter++;
    //         return;
    //     }
    // }

    if (print_info == 1) {
        printf("block_size[id_initial]: %d, np: %ld, id_inital: %d, id: %d, lvl: %d, tree_id: %ld\n", block_size[id_initial], n_points, id_initial, id, lvl, tree_id);
        fflush(stdout);
    }

    // if (lvl == 1) {
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     exit(0);
    // }

    in_min.idx0 = idx_global[0];
    in_min.p_id = id;
    int local = 0;
    for (int i = 1; i < block_size[id_initial]; i++) {
        if (idx_global[i] < in_min.idx0) {
            in_min.idx0 = idx_global[i];
            local = i;
        }
    }

    // printf("id: %d, in_min.idx0: %d, in_min.p_id: %d, local: %d\n", id, in_min.idx0, in_min.p_id, local);
    // fflush(stdout);

    MPI_Allreduce(&in_min, &out_min, 1, MPI_2INT, MPI_MINLOC, comm);

    if (id == out_min.p_id) {
        // printf("idx0: %d ", out_min.p_id);
        // fflush(stdout);
        // print_point(pt_arr[out_min.idx0], n_dims);
        for (int j = 0; j < n_dims; j++) {
            pt_arr_idx0[j] = p_aux[local*n_dims+j];
        }
    }

    MPI_Bcast(pt_arr_idx0, n_dims, MPI_DOUBLE, out_min.p_id, comm);

    // printf("id: %d ", id);
    // fflush(stdout);
    // print_point(pt_arr_idx0, n_dims);

    double d;
    int a = 0;
    in_max.max_d = 0;
    in_max.p_id = id;
    for (int i = 0; i < block_size[id_initial]; ++i) {
        d = 0;
        for (int j = 0; j < n_dims; j++)
            d += (pt_arr_idx0[j] - p_aux[i*n_dims+j]) * (pt_arr_idx0[j] - p_aux[i*n_dims+j]);
        if (d > in_max.max_d) {
            in_max.max_d = d;
            a = i;
        }
    }

    MPI_Allreduce(&in_max, &out_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
    
    if (id == out_max.p_id) {
        for (int j = 0; j < n_dims; j++) {
            pt_arr_a[j] = p_aux[a*n_dims+j];
        }
        // printf("a - id: %d ", id);
        // fflush(stdout);
        // print_point(pt_arr_a, n_dims);
    }

    MPI_Bcast(pt_arr_a, n_dims, MPI_DOUBLE, out_max.p_id, comm);

    int b = 0;
    in_max.max_d = 0;
    in_max.p_id = id;
    for (int i = 0; i < block_size[id_initial]; ++i) {
        d = 0;
        for (int j = 0; j < n_dims; j++)
            d += (pt_arr_a[j] - p_aux[i*n_dims+j]) * (pt_arr_a[j] - p_aux[i*n_dims+j]);
        if (d > in_max.max_d) {
            in_max.max_d = d;
            b = i;
        }
    }

    MPI_Allreduce(&in_max, &out_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);

    if (id == out_max.p_id) {
        for (int j = 0; j < n_dims; j++) {
            pt_arr_b[j] = p_aux[b*n_dims+j];
        }
        // printf("b - id: %d ", id);
        // fflush(stdout);
        // print_point(pt_arr_b, n_dims);
    }

    MPI_Bcast(pt_arr_b, n_dims, MPI_DOUBLE, out_max.p_id, comm);

    if (pt_arr_a[0] > pt_arr_b[0]) {
        double temp;
        for (int i = 0; i < n_dims; i++) {
            temp = pt_arr_a[i];
            pt_arr_a[i] = pt_arr_b[i];
            pt_arr_b[i] = temp;
        }
    }

    if (print_info == 1) {
        printf("id: %d a: ", id);
        fflush(stdout);
        print_point(pt_arr_a, n_dims);

        printf("id: %d b: ", id);
        fflush(stdout);
        print_point(pt_arr_b, n_dims);
    }

    for (int i = 0; i < block_size[id_initial]; ++i) {
        proj[i] = 0;
        for (int j = 0; j < n_dims; j++)
            proj[i] += (p_aux[i*n_dims+j] - pt_arr_a[j]) * (pt_arr_b[j] - pt_arr_a[j]);    
    }

    // MPI_Barrier(MPI_COMM_WORLD);
    // for (int i = 0; i < block_size[id_initial]; i++) {
    //     printf("id: %d, i: %d, idx_global[i]: %ld, proj[i]: %f, pt_arr[i][0]: %f, pt_arr[i][1]: %f\n", id, i, idx_global[i], proj[i], pt_arr[i][0], pt_arr[i][1]);
    //     fflush(stdout);
    // }
    // MPI_Barrier(MPI_COMM_WORLD);

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

    // MPI_Barrier(MPI_COMM_WORLD);
    // printf("id: %d ", id);
    // print_point(pivots, p);
    // MPI_Barrier(MPI_COMM_WORLD);

    MPI_Gather(pivots, p, MPI_DOUBLE, glob_pivots, p, MPI_DOUBLE, 0, comm);

    // MPI_Barrier(MPI_COMM_WORLD);
    // if (!id) {
    //     printf("id: %d ", id);
    //     fflush(stdout);
    //     print_point(glob_pivots, p*p);        
    // }
    // MPI_Barrier(MPI_COMM_WORLD);

    if(!id){
        prev_piv_idx = -1;
        for(int i = 1; i < p; ++i) {
            piv_idx = i * p - prev_piv_idx - 1;
            pivots[i-1] = quickselect_seq_2(piv_idx, prev_piv_idx+1, p*p);
            prev_piv_idx += piv_idx + 1;
        }
    }

    MPI_Bcast(pivots, p-1, MPI_DOUBLE, 0, comm);

    // MPI_Barrier(MPI_COMM_WORLD);
    // printf("id: %d ", id);
    // print_point(pivots, p-1);
    // MPI_Barrier(MPI_COMM_WORLD);

    send_counts[0] = partition(pivots[0], 0, block_size[id_initial]);
    int sum = send_counts[0];
    for(int i = 1; i < p-1; i++) {
        send_counts[i] = partition(pivots[i], sum, block_size[id_initial]) - sum;
        // MPI_Barrier(MPI_COMM_WORLD);
        // printf("id: %d, i: %d, send_counts[i]: %d\n", id, i, send_counts[i]);
        // MPI_Barrier(MPI_COMM_WORLD);
        sum += send_counts[i];
    }
    send_counts[p-1] = block_size[id_initial] - sum;

    // MPI_Barrier(MPI_COMM_WORLD);
    // printf("id: %d ", id);
    // print_point(proj, block_size[id]);
    // MPI_Barrier(MPI_COMM_WORLD);

    // MPI_Barrier(MPI_COMM_WORLD);
    // printf("id: %d, send_counts[0]: %d, send_counts[1]: %d, send_counts[2]: %d, send_counts[3]: %d\n", id, send_counts[0], send_counts[1], send_counts[2], send_counts[3]);
    // MPI_Barrier(MPI_COMM_WORLD);

    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, comm);

    int sum1 = 0;
    int sum2 = 0;
    for(int i = 0; i < p; i++) {
        send_displs[i] = sum2;
        recv_displs[i] = sum1;
        sum1 += recv_counts[i];
        sum2 += send_counts[i];
    }

    if (print_info == 1) {
        printf("id: %d, sum1: %d\n", id, sum1);
        fflush(stdout);
    }

    MPI_Alltoallv(proj, send_counts, send_displs, MPI_DOUBLE, recv_buffer, recv_counts, recv_displs, MPI_DOUBLE, comm);

    double * temp = proj;
    proj = recv_buffer;
    recv_buffer = temp;

    MPI_Alltoallv(idx_global, send_counts, send_displs, MPI_LONG, recv_buffer_long, recv_counts, recv_displs, MPI_LONG, comm);

    long * temp_long = idx_global;
    idx_global = recv_buffer_long;
    recv_buffer_long = temp_long;
    
    for(int i = 0; i < p; ++i){
        send_counts[i] *= n_dims;
        send_displs[i] *= n_dims;
        recv_counts[i] *= n_dims;
        recv_displs[i] *= n_dims;
    }

    MPI_Alltoallv(p_aux, send_counts, send_displs, MPI_DOUBLE, p_aux_2, recv_counts, recv_displs, MPI_DOUBLE, comm);

    // for (int i = 0; i < sum1*n_dims; i++) {
    //     // double t = p_aux_2[i];
    //     printf("p_aux_2[i]: %f\n", p_aux_2[i]);
    //     fflush(stdout);
    // }

    // MPI_Barrier(MPI_COMM_WORLD);
    // exit(0);

    temp = p_aux;
    p_aux = p_aux_2;
    p_aux_2 = temp;

    // MPI_Barrier(MPI_COMM_WORLD);
    // exit(0);

    // for (long i = 0; i < aux; i++)
    //     pt_arr[i] = &p_aux[i * n_dims];

    int size = (id < p/2 + (p%2)) ? sum1 : 0;
    int total_size = 0;

    MPI_Allreduce(&size, &total_size, 1, MPI_INT, MPI_SUM, comm);

    // MPI_Barrier(MPI_COMM_WORLD);
    // printf("id: %d ", id);
    // print_point(proj, sum1);
    // MPI_Barrier(MPI_COMM_WORLD);

    if (print_info == 1) {
        printf("id: %d, size: %d, total_size: %d\n", id, size, total_size);
        fflush(stdout);
    }

    // MPI_Barrier(MPI_COMM_WORLD);
    // if(id == 2) {
    //     for (int i = 0; i < sum1; i++) {
    //         printf("id: %d, i: %d, pt_arr[i][0]: %f, pt_arr[i][1]: %f\n", id, i, pt_arr[i][0], pt_arr[i][1]);
    //         fflush(stdout);
    //     }        
    // }
    // MPI_Barrier(MPI_COMM_WORLD);

    // MPI_Barrier(comm);
    // double mx_proj = proj[0];
    // double m_proj = proj[0];
    // for (int i = 0; i < sum1; i++) {
    //     if (proj[i] > mx_proj)
    //         mx_proj = proj[i];
    //     if (proj[i] < m_proj)
    //         m_proj = proj[i];
    // }
    // printf("id: %d, max_proj: %f, min_proj: %f\n", id, mx_proj, m_proj);
    // fflush(stdout);
    // MPI_Barrier(comm);

    int m1 = -1, m2 = -1;
    double min_proj;
    double u;
    if ((id == p/2) && (p%2 == 1)) {
        m1 = quickselect_seq(n_points/2-!(n_points%2)-(total_size-sum1), 0, sum1);
        if (n_points%2 == 0) {
            min_proj = proj[m1+1];
            m2 = m1+1;
            for (int i = m1+1; i < sum1; i++)
                if (proj[i] < min_proj) {
                    min_proj = proj[i];
                    m2 = i;
                }
            swap(m2, m1+1);
            u = (proj[m1] + proj[m1+1])/2;
            MPI_Send(proj+m1+1, sum1-(m1+1), MPI_DOUBLE, id+1, 0, comm);
            MPI_Send(p_aux+(m1+1)*n_dims, (sum1-(m1+1))*n_dims, MPI_DOUBLE, id+1, 1, comm);
            MPI_Send(idx_global+m1+1, sum1-(m1+1), MPI_LONG, id+1, 2, comm);
            sum1 = m1+1;
        }
        else {
            u = proj[m1];
            MPI_Send(proj+m1, sum1-m1, MPI_DOUBLE, id+1, 0, comm);
            MPI_Send(p_aux+m1*n_dims, (sum1-m1)*n_dims, MPI_DOUBLE, id+1, 1, comm);
            MPI_Send(idx_global+m1, sum1-m1, MPI_LONG, id+1, 2, comm);
            sum1 = m1;
        }
    }
    else if ((id == p/2+1) && (p%2 == 1)) {
        u = 0;
        MPI_Recv(proj+sum1, total_size-(n_points/2), MPI_DOUBLE, id-1, 0, comm, MPI_STATUS_IGNORE);        
        MPI_Recv(p_aux+sum1*n_dims, (total_size-(n_points/2))*n_dims, MPI_DOUBLE, id-1, 1, comm, MPI_STATUS_IGNORE);
        MPI_Recv(idx_global+sum1, total_size-(n_points/2), MPI_LONG, id-1, 2, comm, MPI_STATUS_IGNORE);
        sum1 += total_size-(n_points/2);
    }
    else if ((n_points%2 == 1) && (p%2 == 0)) {
        if (id == p/2 - 1) {
            if (total_size > n_points/2) {
                m1 = quickselect_seq((n_points)/2-(total_size-sum1), 0, sum1);
                u = proj[m1];
                MPI_Send(proj+m1, sum1-m1, MPI_DOUBLE, id+1, 0, comm);
                MPI_Send(p_aux+m1*n_dims, (sum1-m1)*n_dims, MPI_DOUBLE, id+1, 1, comm);
                MPI_Send(idx_global+m1, sum1-m1, MPI_LONG, id+1, 2, comm);
                sum1 = m1;
            }
            else {
                u = 0;
                MPI_Recv(proj+sum1, (n_points/2)-total_size, MPI_DOUBLE, id+1, 0, comm, MPI_STATUS_IGNORE);        
                MPI_Recv(p_aux+sum1*n_dims, ((n_points/2)-total_size)*n_dims, MPI_DOUBLE, id+1, 1, comm, MPI_STATUS_IGNORE);
                MPI_Recv(idx_global+sum1, (n_points/2)-total_size, MPI_LONG, id+1, 2, comm, MPI_STATUS_IGNORE);        
                sum1 += (n_points/2)-total_size;
            }
        }
        else if (id == p/2) {
            if (total_size <= n_points/2) {
                m1 = quickselect_seq((n_points)/2-(total_size), 0, sum1);
                u = proj[m1];
                MPI_Send(proj, m1, MPI_DOUBLE, id-1, 0, comm);
                MPI_Send(p_aux, m1*n_dims, MPI_DOUBLE, id-1, 1, comm);
                MPI_Send(idx_global, m1, MPI_LONG, id-1, 2, comm);
                sum1 = sum1-m1;
                for (int i = 0; i < sum1; i++) {
                    proj[i] = proj[i+m1];
                    idx_global[i] = idx_global[i+m1];
                    for (int j = 0; j < n_dims; j++) {
                        p_aux[i*n_dims+j] = p_aux[(i+m1)*n_dims+j];
                    }
                }
            }
            else {
                u = 0;
                MPI_Recv(proj+sum1, total_size-(n_points/2), MPI_DOUBLE, id-1, 0, comm, MPI_STATUS_IGNORE);        
                MPI_Recv(p_aux+sum1*n_dims, (total_size-(n_points/2))*n_dims, MPI_DOUBLE, id-1, 1, comm, MPI_STATUS_IGNORE);
                MPI_Recv(idx_global+sum1, total_size-(n_points/2), MPI_LONG, id-1, 2, comm, MPI_STATUS_IGNORE);        
                sum1 += total_size-(n_points/2);
            }
        }
        else {
            u = 0;
        }
    }
    else if ((n_points%2 == 0) && (p%2 == 0)) {
        if (id == p/2 - 1) {
            if (total_size >= n_points/2) {
                m1 = quickselect_seq((n_points)/2-1-(total_size-sum1), 0, sum1);
                if (total_size != n_points/2) {
                    min_proj = proj[m1+1];
                    m2 = m1+1;
                    for (int i = m1+1; i < sum1; i++)
                        if (proj[i] < min_proj) {
                            min_proj = proj[i];
                            m2 = i;
                        }
                    swap(m2, m1+1);
                    u = (proj[m1]+proj[m1+1])/2;
                    MPI_Send(proj+m1+1, sum1-(m1+1), MPI_DOUBLE, id+1, 0, comm);
                    MPI_Send(p_aux+(m1+1)*n_dims, (sum1-(m1+1))*n_dims, MPI_DOUBLE, id+1, 1, comm);
                    MPI_Send(idx_global+m1+1, sum1-(m1+1), MPI_LONG, id+1, 2, comm);
                    sum1 = m1+1;
                }
                else {
                    u = proj[m1]/2;
                }
            }
            else {
                u = 0;
                MPI_Recv(proj+sum1, (n_points/2)-total_size, MPI_DOUBLE, id+1, 0, comm, MPI_STATUS_IGNORE);        
                MPI_Recv(p_aux+sum1*n_dims, ((n_points/2)-total_size)*n_dims, MPI_DOUBLE, id+1, 1, comm, MPI_STATUS_IGNORE);
                MPI_Recv(idx_global+sum1, (n_points/2)-total_size, MPI_LONG, id+1, 2, comm, MPI_STATUS_IGNORE);        
                sum1 += (n_points/2)-total_size;
            }
        }
        else if (id == p/2) {
            if (total_size < n_points/2) {
                m1 = quickselect_seq(n_points/2-1-total_size, 0, sum1);
                min_proj = proj[m1+1];
                m2 = m1+1;
                // printf("HERE0 - m1: %d, m2: %d, proj[m1]: %f, proj[m2]: %f, min_proj: %f\n", m1, m2, proj[m1], proj[m2], min_proj);
                // fflush(stdout);
                for (int i = m1+2; i < sum1; i++){
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
                swap(m2, m1+1);
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
                u = (proj[m1]+proj[m1+1])/2;
                MPI_Send(proj, m1+1, MPI_DOUBLE, id-1, 0, comm);
                MPI_Send(p_aux, (m1+1)*n_dims, MPI_DOUBLE, id-1, 1, comm);
                MPI_Send(idx_global, m1+1, MPI_LONG, id-1, 2, comm);
                sum1 = sum1-(m1+1);
                for (int i = 0; i < sum1; i++) {
                    proj[i] = proj[i+m1+1];
                    idx_global[i] = idx_global[i+m1+1];
                    for (int j = 0; j < n_dims; j++) {
                        p_aux[i*n_dims+j] = p_aux[(i+m1+1)*n_dims+j];
                    }
                }
            }
            else if (total_size == n_points/2) {
                m1 = -1;
                min_proj = proj[m1+1];
                m2 = m1+1;
                for (int i = m1+1; i < sum1; i++)
                    if (proj[i] < min_proj) {
                        min_proj = proj[i];
                        m2 = i;
                    }
                swap(m2, m1+1);
                u = proj[0]/2;
            }
            else {
                u = 0;
                MPI_Recv(proj+sum1, total_size-(n_points/2), MPI_DOUBLE, id-1, 0, comm, MPI_STATUS_IGNORE);        
                MPI_Recv(p_aux+sum1*n_dims, (total_size-(n_points/2))*n_dims, MPI_DOUBLE, id-1, 1, comm, MPI_STATUS_IGNORE);
                MPI_Recv(idx_global+sum1, total_size-(n_points/2), MPI_LONG, id-1, 2, comm, MPI_STATUS_IGNORE);
                sum1 += total_size-(n_points/2);
            }
        }
    }
    else {
        u = 0;
    }

    block_size[id_initial] = sum1;

    if (print_info == 1) {
        printf("id: %d, sum1: %d\n", id, sum1);
        fflush(stdout);
    }

    // if (id == 2)
    //     for (int i = 0; i < block_size[id_initial]; i++) {
    //         printf("id: %d, i: %d, idx_global[i]: %ld, proj[i]: %f, pt_arr[i][0]: %f, pt_arr[i][1]: %f\n", id, i, idx_global[i], proj[i], pt_arr[i][0], pt_arr[i][1]);
    //         fflush(stdout);
    //     }

    if (print_info == 1) {
        MPI_Barrier(comm);
        printf("id: %d, u: %f\n", id, u);
        fflush(stdout);
        MPI_Barrier(comm);
    }

    double u_aux;
    MPI_Allreduce(&u, &u_aux, 1, MPI_DOUBLE, MPI_SUM, comm);
    u = u_aux;

    if (print_info == 1) {
        MPI_Barrier(comm);
        printf("id: %d, u: %f\n", id, u);
        fflush(stdout);
        MPI_Barrier(comm);
    }

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
        d = 0;
        for (int j = 0; j < n_dims; j++)
            d += (center[j] - p_aux[i*n_dims+j]) * (center[j] - p_aux[i*n_dims+j]);
        if (d > max_d) {
            max_d = d;
            // printf("id: %d, idx_global[i], : %ld, max_d: %f, pt_arr[i][0]: %f, pt_arr[i][1]: %f\n", id, idx_global[i], max_d, pt_arr[i][0], pt_arr[i][1]);
            // fflush(stdout);
        }
    }

    double rad = 0;

    // MPI_Barrier(MPI_COMM_WORLD);
    // printf("Before - id: %d, max_d: %f\n", id, max_d);
    // fflush(stdout);
    // MPI_Barrier(MPI_COMM_WORLD);

    MPI_Reduce(&max_d, &rad, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

    // MPI_Barrier(MPI_COMM_WORLD);
    // printf("After - id: %d, max_d: %f, rad: %f\n", id, max_d, rad);
    // fflush(stdout);
    // MPI_Barrier(MPI_COMM_WORLD);

    long left = -1;
    long right = -1;

    // printf("%d\n", (int)(log(p_initial)/log(2)-1));

    if (lvl == log(p_initial)/log(2)-1) {
        branch_size[0] = 2 * block_size[0] - 1;
        for (int i = 1; i < p_initial; i++) {
            branch_size[i] = 2 * block_size[i] - 1 + branch_size[i-1];
        }
        // for (int i = 0; i < p; i++) {
        //     printf("id: %d, i: %d, b_size: %d\n", id, i, branch_size[i]);
        //     fflush(stdout);
        // }
        int first_id = pow(2, lvl+1) - 1;
        if (id_initial < 2) {
            left = first_id;
            right = first_id + branch_size[0];
        }
        else if (id == 0){
            left = first_id + branch_size[id_initial-1];
            right = first_id + branch_size[id_initial];
        }
        else {
            left = first_id + branch_size[id_initial-2];
            right = first_id + branch_size[id_initial-1];
        }
    }
    else {
        left = 2 * tree_id + 1;
        right = 2 * tree_id + 2;    
    }

    if (print_info == 1) {
        printf("id: %d, id_initial: %d, left: %ld, right: %ld\n", id, id_initial, left, right);
        fflush(stdout);
    }

    if (!id) {

        tree[tree_counter].node_id = tree_id;
        tree[tree_counter].center_idx = n_center;
        tree[tree_counter].radius = sqrt(rad);
        tree[tree_counter].left = left;
        tree[tree_counter].right = right;

        for (int i = 0; i < n_dims; ++i) {
            centers[n_center][i] = center[i];
        }

        n_center++;

        // printf("%ld %ld %ld %f ", tree_id, tree[tree_counter].left, tree[tree_counter].right, tree[tree_counter].radius);
        // print_point(centers[tree[tree_counter].center_idx], n_dims);
        // fflush(stdout);

        printf("%f ", tree[tree_counter].radius);
        print_point(centers[tree[tree_counter].center_idx], n_dims);
        fflush(stdout);
                
        tree_counter++;
    }

    // if (lvl == 2) {
    //     MPI_Barrier(comm);
    //     exit(0);
    // }

    if (lvl == log(p_initial)/log(2)-1) {

        pts = (pt *)malloc(sizeof(pt) * block_size[id_initial]);
        for (int i = 0; i < block_size[id_initial]; ++i) {
            // pts[i].pt = (double *) malloc(n_dims * sizeof(double));
            // for (int j = 0; j < n_dims; j++) {
            //     pts[i].pt[j] = p_aux[i*n_dims+j];
            // }
            pts[i].pt = &p_aux[i*n_dims];
            // pts[i].pt = p_aux + i * n_dims;
            pts[i].i = idx_global[i];
        }

        if (!id) {
            tree_idx = left;
            ballAlg_seq(0, block_size[id_initial]);
        }
        else {
            tree_idx = right;
            ballAlg_seq(0, block_size[id_initial]);
        }

    }
    else {

        MPI_Comm new_comm;
        // MPI_Barrier(comm);
        MPI_Comm_split(comm, (id+p%2)/(p/2+p%2), id, &new_comm);

        if (tree_id != 0)
            MPI_Comm_free(&comm);

        MPI_Comm_rank(new_comm, &new_id);
        MPI_Comm_size(new_comm, &p);

        if (print_info == 1) {
            printf("id: %d, new_id: %d, tree_id: %ld, proc: %d\n", id, new_id, tree_id, p);
        }

        if (id < p) {
            id = new_id;
            // if (par_info == 1) {
                // printf("id: %d, tree_id: %ld, proc: %d\n", id, tree_id, p);
                // fflush(stdout);
            // }
            ballAlg(n_points/2, left, lvl + 1, new_comm);
        }
        else {
            id = new_id;
            // if (par_info == 1) {
                // printf("id: %d, tree_id: %ld, proc: %d\n", id, tree_id, p);
                // fflush(stdout);
            // }
            ballAlg(n_points/2+n_points%2, right, lvl + 1, new_comm);
        }
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

    n_dims = atoi(argv[1]);
    np = atol(argv[2]);
    unsigned seed = atoi(argv[3]);
    srandom(seed);

    p_initial = p;
    id_initial = id;

    aux = ceil(np/p) * 2;
    p_aux = (double *)malloc(n_dims * aux * sizeof(double));
    // pt_arr = (double **)malloc(aux * sizeof(double *));
    // for (long i = 0; i < aux; i++)
    //     pt_arr[i] = &p_aux[i * n_dims];

    block_size = (int *) malloc (p * sizeof(int));
    branch_size = (int *) malloc (p * sizeof(int));
    for (int i = 0; i < p; i++) {
        block_size[i] = BLOCK_SIZE(i,p,np);
        // printf("id: %d, i: %d, block_size[i]: %d\n", id, i, block_size[i]);
        // fflush(stdout);
    }

    idx_global = (long *) malloc (aux * sizeof(long));

    int BL = BLOCK_LOW(id, p, np);
    int BH = BLOCK_HIGH(id, p, np);
    // printf("id: %d, BL: %d, BH: %d, BS: %d\n", id, BL, BH, block_size[id_initial]);
    // fflush(stdout);
    for (int i = 0; i < np; i++)
        for (int j = 0; j < n_dims; j++) {
            if ((i >= BL) && (i <= BH)) {
                p_aux[(i-BL) * n_dims + j] = RANGE * ((double)random()) / RAND_MAX;
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
    recv_buffer_long = (long *) malloc (aux * sizeof(long));
    p_aux_2 = (double *)malloc(n_dims * aux * sizeof(double));

    tree = (node *)malloc((ceil(log(p_initial)/log(2))+(2*aux-1)) * sizeof(node));
    double *center_aux = (double *)malloc((ceil(log(p_initial)/log(2))+aux-1) * n_dims * sizeof(double));
    centers = (double **)malloc((ceil(log(p_initial)/log(2))+aux-1) * sizeof(double *));
    for (int i = 0; i < (ceil(log(p_initial)/log(2))+aux-1); ++i)
        centers[i] = &center_aux[i * n_dims];
    center = (double *) malloc (n_dims * sizeof(double));

    ballAlg(np, 0, 0, MPI_COMM_WORLD);

    elapsed_time += MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);
    exit(0);

    free(p_aux);
    free(block_size);
    free(branch_size);
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
    free(recv_buffer_long);
    free(p_aux_2);
    free(tree);
    free(center_aux);
    free(centers);
    free(center);
    free(pts);

    if (!id_initial) {
        fprintf(stderr, "%.1lf\n", elapsed_time);
        //print_tree(tree);
    }

    MPI_Finalize();

}