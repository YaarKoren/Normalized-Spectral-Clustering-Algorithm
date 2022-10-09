#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdlib.h>
#include "spkmeans.h"

void initialize_c_array(double** array, PyObject* obj, int array_length, int inner_length);

static PyObject *cal_by_algo(PyObject *self, PyObject *args) {
	/*
	 This function recieves the input goal, the vectors, N and d
	 runs the requested algorithm according to 'goal', and prints the result
	 */
	double **vectors_array, **result_array;
	double *vectors, *result_inner;
	int goal, N, d, error;
	int i;
	PyObject* vectors_obj, *list;

	list = PyList_New(0);
	    if (!list) {
	    	other_error();
	        return NULL;
	    }

    if(!PyArg_ParseTuple(args, "iOii", &goal ,&vectors_obj, &N, &d)) {
    	other_error();
        return NULL;
    }

	/*allocating space for vectors and centroids */

    vectors = (double*)calloc(N*d, sizeof(double));
    vectors_array = (double**)calloc(N, sizeof(double *));

    result_inner = (double*)calloc(N*N, sizeof(double));
   	result_array = (double**)calloc(N, sizeof(double *));

    if(vectors == NULL || vectors_array == NULL || result_array == NULL || result_inner == NULL){
       	free(vectors);
       	free(vectors_array);
       	free(result_array);
       	free(result_inner);
       	other_error();
       }

	/* create matrices for vectors and centroids */
		for (i = 0 ; i < N ; i++) {
			vectors_array[i] = vectors + i*d;
			result_array[i] = result_inner + i*N;
		}

		initialize_c_array(vectors_array, vectors_obj, N, d);

		switch(goal){
		case WAM:{
			cal_wam(result_array, vectors_array,N,d);
			error = 0;
			break;
		}
		case DDG: {
			error = cal_ddg(result_array, vectors_array,N,d);
			break;
		}
		case LNORM:{
			error = cal_lnorm(result_array, vectors_array,N,d);
			break;
		}
		case JACOBI:{
			/*initialize n array for goal JACOBI*/
			double *eigenvalues;
			eigenvalues = (double*)calloc(N, sizeof(double));
			if(eigenvalues == NULL){
				error = 1;
				break;
			}

			error = cal_jacobi(eigenvalues, result_array, vectors_array, N);
			if(error == 0)
				print_vector(eigenvalues, N);
			break;
		}
		default: {
			error = 1;
			break;
		}
	}

		if(error == 1){
			free(result_inner);
			free(result_array);
			free(vectors);
			free(vectors_array);
			other_error();
		}

		print_matrix(result_array, N, N);

		free(vectors);
		free(vectors_array);
		free(result_inner);
		free(result_array);

		  if(error == 1){
		        other_error();
		  }

		return list;
}

static PyObject *cal_spk_vectors(PyObject *self, PyObject *args) {
	/*
	 This function recives the vectors, N, d and the k input
	 calculates and returns the matching T matrix, and a value for k
	 */
	PyObject* vectors_obj;
	double **vectors_array;
	double *vectors;
	double *t_inner, **t_mat;
	int N, d, k, i, j;
    int error;
    PyObject *list, *num, *item;

    list = PyList_New(0);
    if (!list) {
        return NULL;
    }

    /* Parse arguments */
    if(!PyArg_ParseTuple(args, "Oiii", &vectors_obj, &N, &d, &k)) {
        return NULL;
    }

	vectors = (double*)calloc(N*d, sizeof(double));
	vectors_array = (double**)calloc(N, sizeof(double *));

	t_inner = (double*)calloc(N*N, sizeof(double));
	t_mat = (double**)calloc(N, sizeof(double*));

	if(vectors == NULL || vectors_array == NULL
			|| t_inner == NULL || t_mat == NULL){
    	free(vectors);
    	free(vectors_array);
    	free(t_mat);
    	free(t_inner);
    	other_error();
    }

	/* create matrices for vectors and centroids */
	for (i = 0 ; i < N ; i++) {
		vectors_array[i] = vectors + i*d;
		t_mat[i] = t_inner + i*N;
	}

    /* python obj to double** */
	initialize_c_array(vectors_array, vectors_obj, N, d);

	error = spk_alg_make_T_mat(t_mat, vectors_array, N, d, k);

	if( error == -1 || error == 1){
		free(t_mat);
		free(t_inner);
		other_error();
	} else {
		k = error;
	}

	 if (k == 0) {
	     return list;
	}

        free(vectors);
        free(vectors_array);

	/* translating the double** centroids array to a python list [[]] */
     num = PyFloat_FromDouble(0);

     if(list == NULL) {
 		free(t_mat);
     	free(t_inner);
         other_error();
     }

     for(i = 0; i < N; i++) {
         item = PyList_New(k);
         if(item == NULL) {
         	error = 1;
         	break;
         }

         for(j = 0; j < k; j++) {
             PyList_SET_ITEM(item, j,  PyFloat_FromDouble(t_mat[i][j]));
         }

         PyList_Append(list, item);
     }

     num = PyFloat_FromDouble((double) k);
      if (!num) {
           Py_DECREF(list);
           return NULL;
       }
    PyList_Append(list, num);

    free(t_mat);
    free(t_inner);

    return list;
}


static PyObject *fit(PyObject *self, PyObject *args) {
	/*
	 This function receives the vectors, centroids, N, d, k, max_iter and an epsilon
	 and runs the k_means algorithm.
	 the function prints the output centroids, and returns an empty list.
	 */
	PyObject* vectors_obj, *centroids_obj;
	double **vectors_array,  **centroids_array;
	double *vectors, *centroids;
	int N, d, k, max_iter;
	double eps;
    int error;

    int i, j;
    PyObject *list, *item;

    list = PyList_New(0);
        if (!list) {
            return NULL;
        }

    /* Parse arguments */
    if(!PyArg_ParseTuple(args, "OOiiiid", &vectors_obj, &centroids_obj, &N, &d, &k, &max_iter, &eps)) {
        return NULL;
    }

	/*allocating space for vectors and centroids */

	vectors = (double*)calloc(N*d, sizeof(double));
	vectors_array = (double**)calloc(N, sizeof(double *));

	centroids = (double*)calloc(k*d, sizeof(double));
	centroids_array = (double**)calloc(k, sizeof(double *));

	if(vectors == NULL || vectors_array == NULL || centroids == NULL || centroids_array == NULL){
    	free(vectors);
    	free(vectors_array);
	free(centroids);
	free(centroids_array);
    	other_error();
    }

	/* create matrices for vectors and centroids */
	for (i = 0 ; i < N ; i++) {
		vectors_array[i] = vectors + i*d;
	}

	for (i = 0 ; i < k ; i++) {
	    centroids_array[i] = centroids + i*d;
	}

    /* python obj to double** */

       for (i = 0; i < N; i++) {
    	   list = PyList_GetItem(vectors_obj, i);
    	   	for (j = 0; j < d; j++) {
    	   		item = PyList_GetItem(list, j);
    	   		if (!PyFloat_Check(item)){
    	   			Py_DECREF(list);
    	   			return NULL;
				   }
    	   		vectors_array[i][j] = PyFloat_AsDouble(item);
       }
    }

       for (i = 0; i < k; i++) {
		    list = PyList_GetItem(centroids_obj, i);
     	   	for (j = 0; j < d; j++) {
            	item = PyList_GetItem(list, j);
            	if (!PyFloat_Check(item)){
            		Py_DECREF(list);
            		return NULL;
				}
            	centroids_array[i][j] = PyFloat_AsDouble(item);
        	}
     	}

    error = K_means(vectors_array, centroids_array, N, d, k, max_iter, eps);

    free(vectors);
    free(vectors_array);

    if(error == 1){
    	free(centroids);
    	free(centroids_array);
    	other_error();
    }

    print_matrix(centroids_array, k, d);

    free(centroids);
    free(centroids_array);

    if(error == 1){
        other_error();
    }

    return list;
}

void initialize_c_array(double** array, PyObject* obj, int array_length, int inner_length){
	int i, j;
	PyObject *list, *item;

	item = PyFloat_FromDouble(0);
	list = PyList_New(0);

	for (i = 0; i < array_length; i++) {
	    	   list = PyList_GetItem(obj, i);
	    	   	for (j = 0; j < inner_length; j++) {
	    	   		item = PyList_GetItem(list, j);
	    	   		if (!PyFloat_Check(item)){
	    	   			Py_DECREF(list);
	    	   			other_error();
					   }
	    	   		array[i][j] = PyFloat_AsDouble(item);
	       }
	    }
	Py_DECREF(list);
	Py_DECREF(item);
}

static PyMethodDef KMeansMethods[] = {
    {"fit", (PyCFunction) fit, METH_VARARGS, PyDoc_STR("Function that calculates and prints centroids accorning to kMeans algorithm")},
	{"cal_by_algo", (PyCFunction) cal_by_algo, METH_VARARGS, PyDoc_STR("Function that runs algorithm for the input goal, and prints the result")},
	{"cal_spk_vectors", (PyCFunction) cal_spk_vectors, METH_VARARGS, PyDoc_STR("Function that calculates and returns the matching T matrix, and a value for k")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef Moduledef  = {
    PyModuleDef_HEAD_INIT,
    "spkmeansmodule",
    NULL,
    -1,
    KMeansMethods
};

PyMODINIT_FUNC PyInit_spkmeansmodule(void) {
	PyObject *m;

    m = PyModule_Create(&Moduledef);
    if(!m){
    	return NULL;
    }
    return m;
}
