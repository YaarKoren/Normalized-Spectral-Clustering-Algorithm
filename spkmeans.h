
#ifndef SPKMEANS_H_
#define SPKMEANS_H_

#define EPSILON 0.00001
#define JACOBI_MAX_ITER_DEFAULT 100

enum goals {
    WAM,
    DDG,
	LNORM,
    JACOBI
};

static const char* const goal_names[] = {"wam", "ddg", "lnorm", "jacobi" };

/*
vectors - n X d matrix. n vectors, each one is of d coordinates (called X in the assignment descriptoin)
zero_mat - n X n matrix, initialized with 0 in all entries, to store the weighted adjency matrix (wam or W)
the function cahnges zero_mat, makes it the weighted adjacency matrix, out of the n vectors
matrix X -> matrix W
*/
void cal_wam(double** w_matrics, double** vectors, int N, int d);

/*
vectors - n X d matrix. n vectors, each one is of d coordinates (called X in the assignment descriptoin)
zero_mat - n X n matrix, initialized with 0 in all entries, to store the diagonal degree matrix (ddg or D)
the function cahnges zero_mat, makes it the diagonal degree matrix, out of the n vectors
matrix X -> matrix D
returns 0 if successful, 1 if failed
*/
int cal_ddg(double** ddg_matrics, double** vectors, int n, int d);

/*vectors - n X d matrix. n vectors, each one is of d coordinates (called X in the assignment descriptoin)
zero_mat - n X n matrix, initialized with 0 in all entries, to store the L_norm matrix
L_norm = I - D^(-0.5)WD^(-0.5)
the function cahnges zero_mat, makes it the L_norm matrix, out of the n vectors
matrix X -> matrix L_norm
returns 0 if successful, 1 if failed
*/
int cal_lnorm(double** lnorm_matrics, double** vectors, int n, int d);

/*
the func calculate the n eigenvalues and n eigenvectors of sym_mat, a symatric matrix of size n.
changes eighenvalues and eigenvectors, to store them, each vector in a coloumn of the mat eigenvectors.
returns 0 if successful, 1 if failed
*/
int cal_jacobi(double *eigenvalues, double **eigenvectors, double **sym_mat, int n);

/*
steps 1-5 of the spk algorithm - create T mat and calculate k.
zero_mat - n X k matrix, to store the T matrix.
if the input k is 0, calculates k by eigengap heuristic.
returns the value k (the input or the calculated one) if successful, -1 if failed
*/
int spk_alg_make_T_mat(double **zero_mat, double **vectors_mat, int n, int d, int k);

/*utils functions*/

void cal_ddg_from_wam(double **zero_mat, double **wam, int n);

/*
all matrices of size n X n
the fucn computes L_norm matrix and store it in zero_mat
notice that the func changes ddg!
returns 0 if successful, 1 if failed
*/
int cal_lnorm_from_wam_and_ddg(double **zero_mat, double **wam, double **ddg, int n);

void jacobi_alg(double **A, double **new_A, double **prev_P, double **curr_P, double **V, int n);

/*
mat1 - matrix n x d
mat2 - matrix d X m
zero_mat - matrix n X m
the func changes zero_mat to the result of mat1 and mat2 multiplication
*/
void two_matrix_mul(double **zero_mat, double **mat1, double **mat2, int n, int d, int m);

/*
all matrices of size n
the func changes zero_mat to the result of the multiplication: mat1*mat2*mat3
returns 0 if sucesssful, 1 if failed
*/
int three_matrix_mul(double **zero_mat, double **mat1, double **mat2, double **mat3, int n);

/*
mat - diagonal matrix n X d
the func changes mat to mat in the power of exp
*/
void digonal_matrix_exp(double **mat, int n, int d, double exp);

/*makes a matrix n X n Identity matrix*/
void make_I_mat(double **mat, int n);

/*makes a zero matrix n X n Identity matrix*/
void make_I_mat_from_zero_mat(double **zero_mat, int n);

void cal_I_minus_mat(double **zero_mat, double **mat, int n, int d);

void make_list_of_diag_values(double *val_list, double **mat, int n);

/*
A and mat are square matrices of size n
the func changes A to store the values of mat, i.e A = mat
*/
void store_second_in_first(double **A, double **mat, int n);

/*
sym_mat is a symetrical matrix of size n
A is a square matrix of size n
the func changes A to store the values of sym_mat, i.e A = sym_mat
*/
void store_second_in_first_sym(double **A, double **sym_mat, int n);

void transpose_square_mat(double **trans, double **mat, int n);

/*
eigenvalues are sorted in descending order
calculates k by the Eigengap Heuristic, returns k
*/
int cal_k(double *eigenvalues, int n);


/*
evalues - a list on n eigenvalues
evectors - a matrix n X n, each column is an eigenvalue
the eigenvector in column [i]  and the eigenvalue in place [i] are coresponding.
sorted_evalues and sorted_evectors are empty arrays, to store the sorted values and vectors.
the func stores the values in descending order (from big to small)
and the vectors respectivly, each vector in a column.
*/
int sort_evalus_and_evectors(double *sorted_evalues, double **sorted_evectors, double *evalues,
	double **evectors, int n);


/*
sorted eigenvectors - n X n matirx, each column is a vector.
the vectors are sorted from left to right in descending order (big to small), according to the values of their
corresponding eigenvalues.
the func takes the k leftmost vectors, i.e the largest eigenvectors, and from them form a n X k matrix (T mat),
by renormalizing each of U rows to have unit length.
stores the n X k matrix T in zero_mat.
*/
void sorted_eigenvectors_to_T_mat(double **zero_mat, double **sorted_eigenvectors, int n, int k);

struct Eigenval_and_Eigenvector{
	double val;
	double *vector;

};

struct Eigenval_and_Eigenvector* init_Eigenval_and_Eigenvector(double val, double *vector);

/*compartor to sort in descending order*/
int eigenval_comp(const void *p, const void *q);



/*jacobi algorithm*/

/*
finds the indices (i, j) of the off-diagonal element with the largest absolute value,
in a symetric matrix od size n,
and stores them in pivot.
if there are couple of elements with the same abs value, we choose the highest and leftmost.
*/
void find_pivot(int *pivot, double **sym_mat, int n);

/*
calcuates and returns the value of theta in Jacobi algorithm,
in the matrix mat and the indices i and j
*/
double cal_theta(double **mat, int i, int j);

double cal_t_from_theta(double theta);
double cal_c_from_t(double t);
double cal_s_from_t_and_c(double t, double c);

/*
creates the P matrix of Jacobi alogorithm, and store it in zero_mat
the P matrix is of size n, with the elements c and s
*/
void create_P_matrix(double **zero_mat, int n, int i, int j, double c, double s);

/*
calculates and returns the sum of squares of all off-diagonal elements of sym_mat,
a symetrical matrix of size n.
*/
double off(double **sym_mat, int n);

/*
calculates (P_tranpose)AP, accroding to the formula described in the assignment,
using c and s.
the new matrix is stored at new_A
A is a symetrical matrix of size n, after the change new_A is also a symetrical matrix
*/
void cal_newA_from_A(double **new_A, double **A, int n, int i, int j, double c, double s);

/* sym_mat -  symetrical matrix of size n
returns 0 if the matrix is diagonal (zeros in all off-diagonal elements)
1 if not diagonal
*/
int is_not_diagonal_mat(double **sym_mat, int n);

/*kmeans algorithm*/
int K_means(double **vectors_array, double **centroids_array, int N, int d, int K, int max_iter, double eps);
int get_closest_cluster_index(double *v, double **centroids, int d, int K);
void assign_vector_to_cluster (double* v, double **centroids, double ***clusters, int *cluster_size, int d, int K);
void update_centroid_value(double* centroid ,double ** cluster, int cluster_size, int d);
int update_centroids_and_check_distance_difference(double **centroids, double*** clusters, int* cluster_sizes, int k, int d, double eps);
int get_closest_cluster_index(double *v, double **centroids, int d, int K);
void assign_vector_to_cluster (double* v, double **centroids, double ***clusters, int *cluster_size, int d, int K);

/*calculate Euclidean distance of 2 vectors, v1 and v2, both of dimension d*/
double cal_norm (double *v1, double *v2, int d);


/*print results*/

/*print matrix of n X d size*/
void print_matrix(double **matrix, int n, int d);

/*print vector of d dimension*/
void print_vector(double *vector, int d);

/* throw errors */
void other_error();
void input_error();

/*read from file*/

/*get number of rows (vectors)*/
int get_n(FILE *ifp);

/*get number of columns (dimension of vectors)*/
int get_d(FILE *ifp);

int read_from_file(FILE *ifp, double **vectors, int n, int d);

#endif /* SPKMEANS_H_ */

