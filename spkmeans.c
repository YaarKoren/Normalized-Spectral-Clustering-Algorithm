#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "spkmeans.h"

int main(int argc, char **argv){
    int i, n, d, error;
	char *input_path;
    char *goal;
	double *vectors_vals;
	double **vectors_mat;
    double *zero_mat_values;
	double **zero_mat;
	double *eigenvalues;
	FILE *ifp;

    /*read CMD arguments*/
    if(argc != 3){ 
        input_error();
    }

    goal = argv[1];
    input_path = argv[2];


	/*check if goal is legal*/
	if (! (strcmp(goal, goal_names[WAM]) == 0 || strcmp(goal, goal_names[DDG]) == 0 ||
		strcmp(goal, goal_names[LNORM]) == 0 || strcmp(goal, goal_names[JACOBI]) == 0) ){
			input_error();
		}

    /*open file*/
	ifp = fopen(input_path, "r");
	if (ifp == NULL){  
		fclose(ifp);
		input_error();
	}

    n = get_n(ifp);
	rewind(ifp);

    if (strcmp(goal, goal_names[JACOBI]) == 0) { /*strcmp returns 0 iff the strings are equal*/
        d = n;
    }
    else {
        d = get_d(ifp);
	    rewind(ifp);
    }

    /*create an n*d matrix to store the vectors*/
	vectors_vals = (double*)calloc(n*d, sizeof(double));
	vectors_mat = (double**)calloc(n, sizeof(double *));

	if(vectors_vals == NULL || vectors_mat == NULL){
    	free(vectors_vals);
    	free(vectors_mat);
    	other_error();
    }

	for (i = 0 ; i < n ; i++) {
		vectors_mat[i] = vectors_vals + i*d;
	}

    /*read from file and store the vectors*/
	error = read_from_file(ifp, vectors_mat, n, d);
	fclose(ifp);
	if(error) {
		free(vectors_vals);
		free(vectors_mat);
		other_error();
	}

	/*initialize n array for goal JACOBI*/
	eigenvalues = (double*)calloc(n, sizeof(double));

    /*initialize n X n matrix, 0 in all entries*/
    zero_mat_values = (double*)calloc(n*n, sizeof(double));

    /*initialize n-size pointer array, 0 in all entries*/
    zero_mat = (double**)calloc(n, sizeof(double*));

    if(zero_mat_values == NULL || zero_mat == NULL || eigenvalues == NULL){
      	free(zero_mat_values);
    	free(zero_mat);
		free(vectors_vals);
		free(vectors_mat);
		free(eigenvalues);
    	other_error();
    }

    /*connect zero_mat_values and zero_mat, make easy access to zero_mat as a matrix*/
    for (i = 0 ; i < n ; i++) {
    	zero_mat[i] = zero_mat_values + i*n;
	}
	

	if (strcmp(goal, goal_names[WAM]) == 0){
		cal_wam(zero_mat, vectors_mat, n, d);
	}

	else if (strcmp(goal, goal_names[DDG]) == 0){
		error = cal_ddg(zero_mat, vectors_mat, n, d);
	}

	else if (strcmp(goal, goal_names[LNORM]) == 0){
		error = cal_lnorm(zero_mat, vectors_mat, n, d);
	}
	
	else if (strcmp(goal, goal_names[JACOBI]) == 0){
		error = cal_jacobi(eigenvalues, zero_mat, vectors_mat, n);
	}

	if (error) {
		free(zero_mat_values);
    	free(zero_mat);
		free(vectors_vals);
		free(vectors_mat);
		free(eigenvalues);
		other_error();
	}

	if (strcmp(goal, goal_names[JACOBI]) == 0){
		print_vector(eigenvalues, n);
	}

	print_matrix(zero_mat, n, n);

	free(zero_mat_values);
    free(zero_mat);
	free(vectors_vals);
	free(vectors_mat);
	free(eigenvalues);
	exit(0);

}

int cal_jacobi(double *eigenvalues, double **eigenvectors, double **sym_mat, int n){

	double **A;
	double *A_values;
	double **new_A;
	double *new_A_values;
	double **curr_P;
	double *curr_P_values;
	double **prev_P;
	double *prev_P_values;
	int i;

		if (n == 1){
			make_list_of_diag_values(eigenvalues, sym_mat, n);
			eigenvectors[0][0] = 1;
			return 0;
		}

	/*allocate space and create the 4 matrices*/

	/*initialize n X n matrix, 0 in all entries*/
    A_values = (double*)calloc(n*n, sizeof(double));
	new_A_values = (double*)calloc(n*n, sizeof(double));
	curr_P_values = (double*)calloc(n*n, sizeof(double));
	prev_P_values = (double*)calloc(n*n, sizeof(double));

    /*initialize n-size pointer array, 0 in all entries*/
    A = (double**)calloc(n, sizeof(double*));
	new_A = (double**)calloc(n, sizeof(double*));
	curr_P = (double**)calloc(n, sizeof(double*));
	prev_P = (double**)calloc(n, sizeof(double*));

    if(A == NULL || new_A == NULL || curr_P == NULL || prev_P == NULL || A_values == NULL
		 || new_A_values == NULL || curr_P_values == NULL || prev_P_values == NULL ){
      	free(A);
    	free(new_A);
		free(curr_P);
		free(prev_P);
		free(A_values);
		free(new_A_values);
		free(curr_P_values);
		free(prev_P_values);
    	return 1;
    }

    /*connect weights and wam, make easy access to wam as a matrix*/
    for (i = 0 ; i < n ; i++) {
    	A[i] = A_values + i*n;
		new_A[i] = new_A_values + i*n;
		curr_P[i] = curr_P_values + i*n;
		prev_P[i] = prev_P_values + i*n;
	}

	/*store in A copy of sym_mat*/
	store_second_in_first_sym(A, sym_mat, n);

	jacobi_alg(A, new_A, prev_P, curr_P, eigenvectors, n);

	make_list_of_diag_values(eigenvalues, A, n);

	free(A);
    free(new_A);
	free(curr_P);
	free(prev_P);
	free(A_values);
	free(new_A_values);
	free(curr_P_values);
	free(prev_P_values);

	return 0;

}

void jacobi_alg(double **A, double **new_A, double **prev_P, double **curr_P, double **V,
					int n){

	double mat_gap;
	int iter_counter = 0;
	int pivot[2] = {0};
	int i;
	int j;
	double theta;
	double t;
	double c;
	double s;

	/*initialize values before loop*/
	make_I_mat_from_zero_mat(prev_P, n);
	mat_gap = 1;

	while ( (mat_gap > EPSILON) && (iter_counter < JACOBI_MAX_ITER_DEFAULT) ) {

		if ( ! (is_not_diagonal_mat(A , n)) ){ /*A is diagonal*/
			/*eigenvalues are the values on the diagonal*/
			make_I_mat(V, n); /*eigenvectors are the columns of I matrix*/
			return;
		}

		find_pivot(pivot, A, n);
		i = pivot[0];
		j = pivot[1];
		theta = cal_theta(A, i, j);
		t = cal_t_from_theta(theta);
		c = cal_c_from_t(t);
		s = cal_s_from_t_and_c(t, c);

		create_P_matrix(curr_P, n, i, j, c, s);

		/*cal V = V * curr_P*/
		two_matrix_mul(V, prev_P, curr_P, n, n, n);

		/*store V in prev_P
		i.e placement: prev_P = V*/
		store_second_in_first(prev_P, V, n);

		/*cal new_A*/
		cal_newA_from_A(new_A, A, n, i, j, c, s);

		mat_gap = off(A, n) - off(new_A, n);

		/*store new_A in A
		i.e placement: A = new_A*/
		store_second_in_first_sym(A,  new_A, n);

		iter_counter++;

	}

	return;

}

void find_pivot(int *pivot, double **sym_mat, int n){

	int i, j;
	double max_val = -1;

	/*search only in upper-half, so i < j*/
	for (j = 1; j < n; j ++){
		for (i = 0; i < j; i++){
			if (fabs(sym_mat[i][j]) > max_val){
				pivot[0] = i;
				pivot[1] = j;
				max_val = fabs(sym_mat[i][j]);
			}
		}
	}

	return;
}

double cal_theta(double **mat, int i, int j){

	double theta;

	theta = (mat[j][j] - mat[i][i]) / (2 * mat[i][j]);

	return theta;

}

double cal_t_from_theta(double theta){

	double t;
	int sign;

	if (theta == 0 || theta > 0){
		sign = 1;
	}
	else{
		sign = -1;
	}

	t = sign / ( (fabs(theta)) + (pow( (pow(theta, 2) + 1), 0.5)) );

	return t;
}

double cal_c_from_t(double t){

	double c;
	c = 1 / (pow(( pow(t, 2) + 1 ), 0.5));

	return c;
}

double cal_s_from_t_and_c(double t, double c){

	double s;
	s = t * c;
	return s;
}

void create_P_matrix(double **zero_mat, int n, int i, int j, double c, double s){

	make_I_mat(zero_mat, n);
	zero_mat[i][i] = c;
	zero_mat[j][j] = c;
	zero_mat[i][j] = s;
	zero_mat[j][i] = -s;

	return;
}

double off(double **sym_mat, int n){

	int i, j;
	double sum;
	double off;

	sum = 0;
	for (i = 1; i < n; i++){
		for (j = 0; j < i; j++){
			sum += 2* pow(sym_mat[i][j], 2);
		}
	}
	off = sum * 2; /*add twice because the matrix is symetrical*/
	return off;
}

void cal_newA_from_A(double **new_A, double **A, int n, int i, int j, double c, double s){

	int r;
	double val;

	/*initialize new_A to be A
	i.e new_A = A*/
	store_second_in_first(new_A, A, n);


	/*change i and  j rows and columns*/

	/*i row and column*/
	for (r = 0; r < n; r++){
		if ( (r != i) && (r != j ) ){
			val = (c * A[r][i]) - (s * A[r][j]);
			new_A[i][r] = val;
			new_A[r][i] = val;
		}
	}

	/*j row and column*/
	for (r = 0; r < n; r++){
		if ( (r != i) && (r != j ) ){
			val = (c * A[r][j]) + (s * A[r][i]);
			new_A[j][r] = val;
			new_A[r][j] = val;
		}
	}

	new_A[i][i] = ((pow(c, 2)) * A[i][i]) + ((pow(s, 2)) * A[j][j]) - (2 * s * c * A[i][j]);

	new_A[j][j] = ((pow(s, 2)) * A[i][i]) + ((pow(c, 2)) * A[j][j]) + (2 * s * c * A[i][j]);

	new_A[i][j] = 0;
	new_A[j][i] = 0;

	return;
}



void make_list_of_diag_values(double *val_list, double **mat, int n){

	int i;

	for (i = 0; i < n; i++){
		val_list[i] = mat[i][i];
	}

	return;
}


void store_second_in_first_sym(double **A, double **sym_mat, int n){

	int i, j;
	double val;

	/*diagonal*/
	for (i = 0; i < n; i++){
		A[i][i] = sym_mat[i][i];
	}

	/*off-diagonal*/
	for (i = 1; i < n; i++){
		for (j = 0; j < i; j ++){
			val = sym_mat[i][j];
			A[i][j] = val;
			A[j][i] = val;
		}
	}

	return;
}

void store_second_in_first(double **A, double **mat, int n){

	int i, j;

	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			A[i][j] = mat[i][j];
		}
	}

	return;

}

void other_error() {
	printf("An Error Has Occurred\n");
	exit(1);
}

void input_error(){
	printf("Invalid Input!\n");
	exit(1);

}

void print_matrix(double **matrix, int n, int d)
{
    int i;

    for (i = 0; i < n; i++){
		print_vector(matrix[i], d);
    }

	return;
}

void print_vector(double *vector, int d)
{
	int i;

	for (i = 0; i < d; i++){
		if (i < (d-1)){
                printf("%.4f,", vector[i]);
            }
            if (i == (d-1)){
                printf("%.4f\n", vector[i]);
			}
	}
	return;
}

int get_n(FILE *ifp) {

	int count;
	char c;

	count = 0;
	for (c = fgetc(ifp); c != EOF; c = fgetc(ifp))
		if (c == '\n') /*newline*/
			count++;

	/*last line is always empty, but in case of 4 lines there will be exactly 4 '\n's*/
	return count;
}

int get_d(FILE *ifp) {
	int count;
	char c;

	count = 0;
	while  ( (c = fgetc(ifp) ) != '\n') {
		if (c == ',')
			count++;
	}

	count++; /*because the line ends with "\n" not ","*/
	return count;

}

int read_from_file(FILE *ifp, double **vectors, int n, int d) {
	int i, j;
	char c;
	char *start_of_str_coor, *str_coor;
	double coor;

	/*store the vectors*/
	str_coor = (char *)malloc(20 * sizeof(char)); /* 7 becausse longest num i.e.: "-7.8602"*/
	if (str_coor == NULL){
		free(str_coor);
		return 1;
	}
	for  (i = 0; i < n; i++) { /* run n times - number of vectors = number of lines (minus the empty one)*/
		for (j = 0; j < d; j++) { /* run d times - number of coordinates = number of numbers in a line*/
			start_of_str_coor = str_coor;
			while ( (c=fgetc(ifp)) != ',' && c != '\n') { /*collect the chars of the coordinate*/
				*str_coor = (char)c;
				str_coor++;
			}
			coor = strtod(start_of_str_coor, NULL);
			vectors[i][j] = coor;
			str_coor = start_of_str_coor;
		}

	}
	free(str_coor);
	return 0;
}

double cal_norm (double *v1, double *v2, int d) {

	int i;
	double sum = 0, norm;
	for (i=0; i<d; i++){
		sum += pow ( (*(v1+i)) - (*(v2+i)), 2 );
	}

	norm = sqrt (sum);
	return norm;
}

void cal_wam (double **zero_mat, double **vectors, int n, int d){

	double norm, exponent, value;
    int i, j;

    for (i = 1 ; i < n ; i++){
        for (j = 0; j < i; j++){
            norm = cal_norm(vectors[i], vectors[j], d);
            exponent = (-1) * (norm/2);
			value = exp(exponent); /*the exp func returns e ^ exponent*/
            zero_mat[i][j] = value; /*fill the lower half of wam (below the main diagonal)*/
			zero_mat[j][i] = value; /*fill the upper half of wam (above the main diagonal)*/
			/*the matix is symetrical*/
			/*the main diagonal remains zeros*/
        }
    }
    return;
}

void cal_ddg_from_wam(double **zero_mat, double **wam, int n){

	int i, j;
	double sum;

	for (i = 0; i < n; i++){
		sum = 0;
		for (j = 0; j < n; j++){
			sum += wam[i][j];
		}
		zero_mat[i][i] = sum;
		sum = 0;
	}

	return;
}

int cal_ddg (double **zero_mat, double **vectors, int n, int d){
	
	double **wam;
    double *weights;
	int i;

	/*create matrix to store wam*/
	
	/*initialize n X n matrix, 0 in all entries*/
    weights = (double*)calloc(n*n, sizeof(double));

    /*initialize n-size pointer array, 0 in all entries*/
    wam = (double**)calloc(n, sizeof(double*));

    if(weights == NULL || wam == NULL){
      	free(weights);
    	free(wam);
    	return 1;
    }

    /*connect weights and wam, make easy access to wam as a matrix*/
    for (i = 0 ; i < n ; i++) {
    	wam[i] = weights + i*n;
	}
	
	cal_wam(wam, vectors, n, d);

	cal_ddg_from_wam(zero_mat, wam, n);
	
	free(weights);
    free(wam);

	return 0;
}


int cal_lnorm_from_wam_and_ddg(double **zero_mat, double **wam, double **ddg, int n){

	double **mul;
	double *mul_values;
	int i, error;

	/*make ddg (D) ddg^(-0.5)*/	
	digonal_matrix_exp(ddg, n, n, -0.5);


	/*create matrix to store mul*/
	
	/*initialize n X n matrix, 0 in all entries*/
    mul_values = (double*)calloc(n*n, sizeof(double));

    /*initialize n-size pointer array, 0 in all entries*/
    mul = (double**)calloc(n, sizeof(double*));

    if(mul_values == NULL || mul == NULL){
      	free(mul_values);
    	free(mul);
    	return 1;
    }

    /*connect mul_values and mul, make easy access to mul as a matrix*/
    for (i = 0 ; i < n ; i++) {
    	mul[i] = mul_values + i*n;
	}

	error = three_matrix_mul(mul, ddg, wam, ddg, n);
	if (error){
		free(mul_values);
    	free(mul);
		return 1;
	}

	cal_I_minus_mat(zero_mat, mul, n, n);

	free(mul_values);
    free(mul);

	return 0;
}


int cal_lnorm(double **zero_mat, double **vectors, int n, int d){

	double **wam;
    double *weights;
	double **ddg;
	double *ddg_values;
	int i, error;

	/*calculate W (wam)*/
	
	/*create matrix to store wam*/
	
	/*initialize n X n matrix, 0 in all entries*/
    weights = (double*)calloc(n*n, sizeof(double));

    /*initialize n-size pointer array, 0 in all entries*/
    wam = (double**)calloc(n, sizeof(double*));

    if(weights == NULL || wam == NULL){
      	free(weights);
    	free(wam);
    	return 1;
    }

    /*connect weights and wam, make easy access to wam as a matrix*/
    for (i = 0 ; i < n ; i++) {
    	wam[i] = weights + i*n;
	}
	
	cal_wam(wam, vectors, n, d);

	/*calculte D (ddg)*/

	/*create matrix to store ddg*/
	
	/*initialize n X n matrix, 0 in all entries*/
    ddg_values = (double*)calloc(n*n, sizeof(double));

    /*initialize n-size pointer array, 0 in all entries*/
    ddg = (double**)calloc(n, sizeof(double*));

    if(ddg_values == NULL || ddg == NULL){
      	free(ddg_values);
    	free(ddg);
    	return 1;
    }

    /*connect ddg_values and ddg, make easy access to ddg as a matrix*/
    for (i = 0 ; i < n ; i++) {
    	ddg[i] = ddg_values + i*n;
	}
	
	cal_ddg_from_wam(ddg, wam, n);

	error = cal_lnorm_from_wam_and_ddg(zero_mat, wam, ddg, n);

	free(weights);
    free(wam);
	free(ddg_values);
    free(ddg);

	if (error){
		return 1;
	}

	return 0;
}

void two_matrix_mul(double **zero_mat, double **mat1, double **mat2, int n, int d, int m){

	int i, j, k;
	double sum;

	for (i= 0; i < n; i++){
		for (j = 0; j < m; j++){
			sum = 0;
			for (k = 0; k < d; k++){
				sum += (mat1[i][k] * mat2[k][j]);
			} 
			zero_mat[i][j] = sum;
		}
	}

	return;
}

int three_matrix_mul(double **zero_mat, double **mat1, double **mat2, double **mat3, int n){

	double **mul;
	double *mul_values;
	int i;

	/*create matrix to store mul*/
	
	/*initialize n X n matrix, 0 in all entries*/
    mul_values = (double*)calloc(n*n, sizeof(double));

    /*initialize n-size pointer array, 0 in all entries*/
    mul = (double**)calloc(n, sizeof(double*));

    if(mul_values == NULL || mul == NULL){
      	free(mul_values);
    	free(mul);
    	return 1;
    }

    /*connect mul_values and mul, make easy access to mul as a matrix*/
    for (i = 0 ; i < n ; i++) {
    	mul[i] = mul_values + i*n;
	}

	two_matrix_mul(mul, mat1, mat2, n, n, n);

	two_matrix_mul(zero_mat, mul, mat3, n, n, n);

	free(mul_values);
    free(mul);
    return 0;
}

void digonal_matrix_exp(double **mat, int n, int d, double exp){

	int i;

	for (i = 0; ( (i <  n) & (i < d) ); i++){
		mat[i][i] = pow(mat[i][i], exp);
	}

	return;
}

void make_I_mat(double **mat, int n){

	int i, j;

	for (i= 0; i < n; i++){
		for (j = 0; j < n; j++){
			if (i == j){
				mat[i][j] = 1;
			}
			else{
				mat[i][j] = 0;
			}
		}
	}
	return;
}

void make_I_mat_from_zero_mat(double **zero_mat, int n){

	int i;
	for (i= 0; i < n; i++){
		zero_mat[i][i] = 1;
	}
	return;
}

void cal_I_minus_mat(double **zero_mat, double **mat, int n, int d){

	int i, j;
	
	for (i = 0; i < n; i ++){
		for (j = 0; j < d; j ++){
			if (i == j){
				zero_mat[i][j] = 1 - mat[i][j];
			}
			else{
				zero_mat[i][j] = - mat[i][j];
			}

		}
	}
	return;
}

int K_means(double **vectors_array, double **centroids_array, int N, int d, int K, int max_iter, double eps){

	int i, iteration_number, flag;
	double **clusters;
	double ***clusters_array;
	int *clusters_sizes;

    /*allocate space for clusters. each cluster is an array of pointers
    (that is, clusters_array in the [i][j] contains pointer to a vector)*/
    clusters = (double **)calloc(K, N*(sizeof(double *)));
    clusters_array = (double***)calloc(K, sizeof(double **));
    /*allocate space for clusters sizes*/
    clusters_sizes = (int*)calloc(K, sizeof(int));

    if(clusters == NULL || clusters_array == NULL || clusters_sizes == NULL){
    	free(clusters);
    	free(clusters_array);
    	free(clusters_sizes);
    	return 1;
    }

    for (i = 0 ; i < K ; i++) { /*make K arrays, each for one cluster*/
        clusters_array[i] = clusters + i*N;
    }

    /*repeat until convergence or until iteration number = max iter)*/
    iteration_number = 0;
    while (iteration_number < max_iter) {

		/*reser cluster sizes*/
		for (i = 0 ; i < K ; i++){
    		clusters_sizes[i] = 0;
    	}

		/*for each vector (N in total), assign vector to the closest cluster*/
        for (i = 0 ; i < N ; i++){
            /*the vector we want is vectors_array[i]*/
            assign_vector_to_cluster (vectors_array[i], centroids_array, clusters_array, clusters_sizes, d, K);
        }

        /*for each centroid (K in total), update centroid, and check if the change is < EPSILON*/
        flag = update_centroids_and_check_distance_difference(centroids_array, clusters_array, clusters_sizes, K, d, eps);

        /*if there was an error in update_centroids function*/
        if(flag <0){
        	free(clusters);
        	free(clusters_array);
        	free(clusters_sizes);
        	return 1;
        }

        /*check convergence
        flag == 0 iff all changes are smaller then EPSILON, that is, there is convergence*/
        if (flag == 0){
            break;
        }

        iteration_number++;
    }

    free(clusters);
    free(clusters_array);
    free(clusters_sizes);

    return 0;
 }

int get_closest_cluster_index(double *v, double **centroids, int d, int K) {
	int i;
	int closest_cluster_index = -1;
	double min_distance, distance;
	for (i=0; i<K; i++){
		distance = cal_norm(v, (*(centroids+i)), d); /*the centroid we check is *(centorids+i)*/
		if (closest_cluster_index < 0 || distance < min_distance) {
			min_distance = distance;
			closest_cluster_index = i;
		}
	}
	return closest_cluster_index;
}


void assign_vector_to_cluster (double* v, double **centroids, double ***clusters, int *cluster_size, int d, int K){
	int cluster_index, new_vector_index;

	/*find the closest cluster*/
	cluster_index = get_closest_cluster_index(v, centroids,d, K);

    		/*check how many vectors already in cluster, to add v in the right place
    		num_of_vectors_in_cluster is *(cluster_size+i)
    		if we had x vectors, places 0, 1, ... , x-1 are occupid, and new vector goes to place x*/
    new_vector_index = cluster_size[cluster_index];

	/*assign the vector
	our cluster is = *(clusters + cluster_index)
	the place to put the vector is = *(cluster + new_vector_index)*/
	clusters[cluster_index][new_vector_index] = v;

	/*update cluster size*/
	cluster_size[cluster_index] ++;

}

void update_centroid_value(double* centroid ,double ** cluster, int cluster_size, int d){
	int i,j;

	/*initialize the centroid*/
	for(i=0;i<d;i++) {
		centroid[i] = 0;
	}

	if(cluster_size == 0) {
		return;
	}

	/*calculate new value for centroid*/
	for(i=0;i<cluster_size;i++) {
		for(j=0;j<d;j++) {
			centroid[j] +=cluster[i][j];
		}
	}
	for(i=0;i<d;i++) {
		centroid[i] /= cluster_size;
	}
}



int update_centroids_and_check_distance_difference(double **centroids, double*** clusters, int* cluster_sizes, int k, int d, double eps){
	int flag, i, j;
	double norm;
	double *temp_vector;

	flag = 0;
	temp_vector = (double*)malloc(d*sizeof(double));

	if(temp_vector == NULL){
		free(temp_vector);
		return -1;
	}

	for(i=0;i<k;i++){ /*for every centroid, update*/
		for(j=0;j<d;j++){ /*save old value of centroid*/
			temp_vector[j] = centroids[i][j];
		}
		update_centroid_value(centroids[i], clusters[i], cluster_sizes[i], d);
		norm = cal_norm(temp_vector, centroids[i],d);
		if(norm >= eps){
			flag =1;
		}
	}
	free(temp_vector);
	return flag;
}


void transpose_square_mat(double **trans, double **mat, int n){

	int i, j;

	for(i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			trans[i][j] = mat[j][i];
		}
	}

	return;
}



struct Eigenval_and_Eigenvector* init_Eigenval_and_Eigenvector(double val, double *vector) {

	struct Eigenval_and_Eigenvector* element = (struct Eigenval_and_Eigenvector*)calloc(1, sizeof(struct Eigenval_and_Eigenvector));
	element->val = val;
	element->vector = vector;
  	return element;
}



int eigenval_comp(const void *p, const void *q){

	double a;
	double b;

	a = ((struct Eigenval_and_Eigenvector *)p)->val;
    b = ((struct Eigenval_and_Eigenvector *)q)->val;

	if (a == b){
		return 0;
	}
	else if (a < b){
		return 1;
	}

	return -1;

}

int sort_evalus_and_evectors(double *sorted_evalues, double **sorted_evectors, double *evalues,
	double **evectors, int n){

	int i;
	double *vectors_values;
	double **vectors;
	struct Eigenval_and_Eigenvector *vals;

	/*make a matrix in which, each row is an eigenvector*/
	/*initialize n X n matrix, 0 in all entries*/
    vectors_values = (double*)calloc(n*(n), sizeof(double));

    /*initialize n-size pointer array, 0 in all entries*/
    vectors = (double**)calloc(n, sizeof(double*));

    if(vectors_values == NULL || vectors == NULL){
      	free(vectors_values);
    	free(vectors);
    	return 1;
    }

    /*connect vectors_values and vectors, make easy access to vectors as a matrix*/
    for (i = 0 ; i < n ; i++) {
    	vectors[i] = vectors_values + i*(n);
	}

	transpose_square_mat(vectors, evectors, n);

	/*create an array of val&vectors, to sort them by vals*/

	vals = (struct Eigenval_and_Eigenvector*)calloc(n, sizeof(struct Eigenval_and_Eigenvector));

	for (i = 0; i < n; i++){
		vals[i] = * init_Eigenval_and_Eigenvector(evalues[i], vectors[i]);
	}

	qsort(vals, n, sizeof(struct Eigenval_and_Eigenvector), eigenval_comp);

	/*store sorted elemnts back*/
	/*eigenvalues*/
	for (i = 0 ; i < n ; i++) {
    	sorted_evalues[i] = vals[i].val;
	}

	/*eigenvectors*/
	for (i = 0 ; i < n ; i++) {
    	vectors[i] = vals[i].vector;
	}
	/*make the vectors as columns*/
	transpose_square_mat(sorted_evectors, vectors, n);

	return 0;

}

int spk_alg_make_T_mat(double **zero_mat, double **vectors_mat, int n, int d, int k){

	int i, error;
	int new_k;

	double *lnorm_vals;
	double **lnorm;
	double *eigenvectors_vals;
	double **eigenvectors;
	double *eigenvalues;
	double *sorted_eigenvectors_vals;
	double **sorted_eigenvectors;
	double *sorted_eigenvalues;


	eigenvalues = (double*)calloc(n, sizeof(double));
	sorted_eigenvalues = (double*)calloc(n, sizeof(double));

	lnorm_vals = (double*)calloc(n*n, sizeof(double));
	eigenvectors_vals = (double*)calloc(n*n, sizeof(double));
	sorted_eigenvectors_vals = (double*)calloc(n*n, sizeof(double));

	lnorm = (double**)calloc(n, sizeof(double *));
	eigenvectors = (double**)calloc(n, sizeof(double *));
	sorted_eigenvectors = (double**)calloc(n, sizeof(double *));

	if(lnorm_vals == NULL || lnorm == NULL ||eigenvectors_vals == NULL || eigenvectors == NULL
		 || eigenvalues == NULL || sorted_eigenvalues == NULL || sorted_eigenvectors_vals == NULL
		 || sorted_eigenvectors == NULL){
		free(lnorm_vals);
		free(lnorm);
		free(eigenvectors_vals);
		free(eigenvectors);
		free(eigenvalues);
		free(sorted_eigenvalues);
		free(sorted_eigenvectors_vals);
		free(sorted_eigenvectors);
    	return -1;
    }

	for (i = 0 ; i < n ; i++) {
		lnorm[i] = lnorm_vals + i*n;
		eigenvectors[i] = eigenvectors_vals + i*n;
		sorted_eigenvectors[i] = sorted_eigenvectors_vals  + i*n;
	}

	error = cal_lnorm(lnorm, vectors_mat, n, d);
	if (error){
		free(lnorm_vals);
		free(lnorm);
		free(eigenvectors_vals);
		free(eigenvectors);
		free(eigenvalues);
		free(sorted_eigenvalues);
		free(sorted_eigenvectors_vals);
		free(sorted_eigenvectors);
    	return -1;
    }

	/*find eigenvalues and eigenvectors of lnorm*/
	error = cal_jacobi(eigenvalues, eigenvectors, lnorm, n);
	if (error){
		free(lnorm_vals);
		free(lnorm);
		free(eigenvectors_vals);
		free(eigenvectors);
		free(eigenvalues);
		free(sorted_eigenvalues);
		free(sorted_eigenvectors_vals);
		free(sorted_eigenvectors);
    	return -1;
    }

	/*sort eigenvalues and eigen vectors*/
	error = sort_evalus_and_evectors(sorted_eigenvalues, sorted_eigenvectors, eigenvalues,
		eigenvectors, n);

	if (error){
		free(lnorm_vals);
		free(lnorm);
		free(eigenvectors_vals);
		free(eigenvectors);
		free(eigenvalues);
		free(sorted_eigenvalues);
		free(sorted_eigenvectors_vals);
		free(sorted_eigenvectors);
    	return -1;
    }

	if (k==0){
		new_k = cal_k(sorted_eigenvalues, n);
	}
	else{
		new_k = k;
	}

	sorted_eigenvectors_to_T_mat(zero_mat, sorted_eigenvectors, n, new_k);

	free(lnorm_vals);
	free(lnorm);
	free(eigenvectors_vals);
	free(eigenvectors);
	free(eigenvalues);
	free(sorted_eigenvalues);
	free(sorted_eigenvectors_vals);
	free(sorted_eigenvectors);

	return new_k;
}


int cal_k(double *eigenvalues, int n){
	int i, k;
	double gap, max_gap;
	max_gap = 0;
	k=0;
		for (i = 0; i<n/2; i++){
			gap = fabs(eigenvalues[i] - eigenvalues[i+1]);
			if (gap > max_gap){
				max_gap = gap;
				k = i + 1;
			}
		}
		return k;
}



void sorted_eigenvectors_to_T_mat(double **zero_mat, double **sorted_eigenvectors, int n, int k){

	int i, j, r;
		double denominator;
		double sum;

		for (i = 0; i < n; i ++){
			sum = 0;
			for (j = 0; j < k; j++){
				sum += pow(sorted_eigenvectors[i][j] ,2);
				denominator = pow(sum, 0.5);
			}
			for (r = 0; r < k; r++){
				if (denominator == 0){
					zero_mat[i][r] = 0;
				} else {
					zero_mat[i][r] = sorted_eigenvectors[i][r] / denominator;
				}
			}
		}
		return;

}

int is_not_diagonal_mat(double **sym_mat, int n) {
	int i, j;
	/*check the upper half of the matrix*/
	for (j = 1; j < n; j ++){
		for (i = 0; i < j; i++){
			if (sym_mat[i][j] != 0){
				return 1;
			}
		}
	}
	return 0;
}
