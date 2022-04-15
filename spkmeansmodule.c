#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <stdlib.h>
#include "spkmeans.h"

typedef double* point;

typedef double** matrix;


struct cluster {
    point  sum;
    int    size;
} ;

struct program_args {
    int K;
    int max_iter;
    double eps;
} ;

struct datapoint_info {
    int N;
    int d;
} ;

struct rotation_mat_info {
    int i;  int j;  double c;  double s;
} ;

static void matrix_free(matrix *mat, int num_of_rows);



/*
----------------------------------------------------------------------- 
-----------------------------EX2 CODE---------------------------------- 
-----------------------------------------------------------------------
*/



/* populate an array of points based on a python-object */
/* input: py_obj is a 1d-list, representing a 2d-list with dimensions dim1 x dim2 */
/* read_input is "free"-safe! Even if it returns -1, there's no mess to clean up */
int read_input(PyObject *py_obj, point **c_arr, int dim1, int dim2)
{
    int i,j;
    double tmp;

    /* allocate memory to datapoint array */
    *c_arr = malloc(sizeof(point) * dim1);
    if (c_arr == NULL)
    {
        return -1;
    }
    
    /* read data */
    for (i = 0; i < dim1; i++)
    {
        /* allocate memory to the next datapoint */
        (*c_arr)[i] = malloc(sizeof(double) * dim2);
        if ((*c_arr)[i] == NULL)
        {
            /* free the memory it already allocated, then finish */
            matrix_free(c_arr, i);
            return -1;
        }

        for (j = 0; j < dim2; j++)
        {
            tmp = PyFloat_AsDouble(PyList_GetItem(py_obj, i * dim2 + j));
            (*c_arr)[i][j] = tmp;
        }        
    }
    return 0;
}

void assign_to_cluster(point xi, struct cluster* clusters, point* mu, int K, int d)
{
    int j;
    int argmin = 0;
    double curr_val;
    double min_val = norm_of_diff(xi, mu[0], d, 0);

    /* find right cluster for xi */
    for (j = 1; j < K; j++)
    {
        curr_val = norm_of_diff(xi, mu[j], d, 0);
        if (curr_val < min_val)
        {
            argmin = j;
            min_val = curr_val;
        }
    }
    
    /* add xi to cluster (update "vector_sum" and "size") */
    for (j = 0; j < d; j++)
    {
        clusters[argmin].sum[j] += xi[j];
    }
    clusters[argmin].size++;
}

int not_converge(point* mu, point* prev_mu, double eps, int K, int d)
{
    int i;
    for (i = 0; i < K; i++)
    {
        if (norm_of_diff(mu[i], prev_mu[i], d, 1) >= eps)
        {
            return 1;
        }
    }
    
    return 0;
}

void free_clusters(struct cluster* clusters, int fail_index)
{
    int i;
    for (i = 0; i < fail_index; i++)
    {
        free(clusters[i].sum);
    }
    free(clusters);
}

static int get_clusters(point *datapoints, point *mu, int N, int d, int K, int max_iter, double eps)
{
    /* declare variables */
    int     i,j;
    int     curr_iter = 0;
    point   *prev_mu;
    
    struct datapoint_info info = {N, d};
    struct program_args   argum = {K, max_iter, eps};
    struct cluster        *clusters;
    

    /* --------------------- NON-INPUT MEMORY ALLOCATION ---------------------------- */
    /* ---------------------    (mu is given as input)   ---------------------------- */

    clusters = malloc(argum.K * sizeof(struct cluster));
    if (clusters == NULL)
    {
        printf("An Error Has Occurred\n");
        matrix_free(&datapoints, info.N); /* free_datapoints(datapoints, info.N); */
        matrix_free(&mu, K);
        return 1;
    }

    for (i = 0; i < argum.K; i++)
    {
        clusters[i].sum = malloc(info.d * sizeof(double));
        if (clusters[i].sum == NULL)
        {            
            printf("An Error Has Occurred\n");
            matrix_free(&datapoints, info.N);  /* free_datapoints(datapoints, info.N); */
            matrix_free(&mu, K);
            free_clusters(clusters, i);
            return 1;
        }
    }

    
    prev_mu = malloc(argum.K * sizeof(point));
    if (prev_mu == NULL)
    {            
        printf("An Error Has Occurred\n");
        matrix_free(&datapoints, info.N);  /* free_datapoints(datapoints, info.N); */
        free_clusters(clusters, argum.K);
        matrix_free(&mu, K); /* free_mu(mu, prev_mu, K, -1); */
        return 1;
    }

    for (i = 0; i < argum.K; i++)
    {
        prev_mu[i] = malloc(info.d * sizeof(double));
        if (prev_mu[i] == NULL)
        {            
            printf("An Error Has Occurred\n");
            matrix_free(&datapoints, info.N);  /* free_datapoints(datapoints, info.N); */
            free_clusters(clusters, argum.K);
            matrix_free(&mu, K); /* free_mu(mu, prev_mu, K, i); */
            matrix_free(&prev_mu, i); /* free_datapoints(prev_mu, i); */
            return 1;
        }   
    }

    /* initialize prev_mu */
    for (i = 0; i < argum.K; i++)
    {
        for (j = 0; j < info.d; j++)
        {
            prev_mu[i][j] = 0;
        }
    }
    

    /* --------------------- PERFORM ALGORITHM ---------------------------- */


    /* perform algorithm */
    while (curr_iter < argum.max_iter && not_converge(mu, prev_mu, argum.eps, argum.K, info.d))
    {
       /* reinitialize clusters */
        for (i = 0; i < argum.K; i++)
        {
            clusters[i].size = 0;
            for (j = 0; j < info.d; j++)
            {
                clusters[i].sum[j] = 0;
            }
        }

        /* distribute points to clusters */
        for (i = 0; i < info.N; i++)
        {
            assign_to_cluster(datapoints[i], clusters, mu, argum.K, info.d);
        }

        /* update mu, prev_mu */
        for (i = 0; i < argum.K; i++)
        {
            for (j = 0; j < info.d; j++)
            {
                prev_mu[i][j] = mu[i][j];

                /* BEWARE of DIVISION BY ZERO */
                if ((double)clusters[i].size == 0)
                    mu[i][j] = 0;
                else
                    mu[i][j] = clusters[i].sum[j] / (double)clusters[i].size;
            }
        }

        curr_iter++;
   }

    /* free memory EXCEPT MU (will be freed in "fit") */
    matrix_free(&datapoints, info.N);  /* free_datapoints(datapoints, info.N); */
    free_clusters(clusters, argum.K);
    matrix_free(&prev_mu, argum.K);  /* free_mu(NULL, prev_mu, 0, argum.K); */

    return 0;
}

static PyObject* fit(PyObject *self, PyObject *args)
{
    PyObject *initial_mu_obj;
    PyObject *datapoints_obj;
    point *initial_mu;
    point *datapoints;
    int N,d,K, max_iter;
    double eps;
    int i,j;

    /* parse arguments from python int our variables */
    if (!PyArg_ParseTuple(args, "OOiiiid", &initial_mu_obj, &datapoints_obj, &N, &d, &K, &max_iter, &eps))
    {
        printf("An Error Has Occurred\n");
        exit(1);  
    }

    /* copy input from python objects into c arrays 
       (read_input returns -1 if function failed) */
    if (read_input(initial_mu_obj, &initial_mu, K, d) == -1 || 
        read_input(datapoints_obj, &datapoints, N, d) == -1)
    {
        if (initial_mu != NULL)
            matrix_free(&initial_mu, K);
        if (datapoints != NULL)
            matrix_free(&datapoints, N);

        printf("An Error Has Occurred\n");
        exit(1);
    }


    /* perform algorithm */
    if (get_clusters(datapoints, initial_mu, N, d, K, max_iter, eps) == 1)
        exit(1);


    /* return chosen centroids */
    /* code is based on this stack-overflow thread:
        https://stackoverflow.com/questions/50668981/how-to-return-a-list-of-ints-in-python-c-api-extension-with-pylist/50683462 */
    PyObject* returned_list = PyList_New(K);
    for (i = 0; i < K; ++i)
    {
        PyObject* mu_i = PyList_New(d);
        for (j = 0; j < d; ++j)
        {
            PyObject* python_double = Py_BuildValue("d", initial_mu[i][j]);
            PyList_SetItem(mu_i, j, python_double);
        }
        PyList_SetItem(returned_list, i, mu_i);
    }

    matrix_free(&initial_mu, K);  /* free_mu(initial_mu, NULL, K, 0); */

    return returned_list;
}



/*
----------------------------------------------------------------------- 
-----------------------------NEW CODE---------------------------------- 
-----------------------------------------------------------------------
*/


/*  returns a matrix initialized with zeroes  */
/*  DDG functions assumes the matrix has zeroes in it !!!  */
static int matrix_malloc(matrix *mat, int dim1, int dim2)
{
    int i,j;

    (*mat) = malloc(sizeof(point) * dim1);
    if ((*mat) == NULL)
    {
        return -1;
    }
    for (i = 0; i < dim1; i++)
    {
        (*mat)[i] = malloc(sizeof(double) * dim2);

        /* if allocation failed - free all allocated arrays, then return -1
           else - populate array with zeroes */
        if ((*mat)[i] == NULL)
        {
            for (j = 0; j < i; j++) 
            {
                free((*mat)[j]);
            }
            free((*mat));
            return -1;
        }
        else
        {
            for (j = 0; j < dim2; j++)
            {
                (*mat)[i][j] = 0;
            }
        }
    }
    return 0;
}

static PyObject* matrix_to_python(matrix mat, int dim1, int dim2)
{
    /* code is based on this stack-overflow thread:
        https://stackoverflow.com/questions/50668981/how-to-return-a-list-of-ints-in-python-c-api-extension-with-pylist/50683462 */
        
    int i,j;
    PyObject* returned_list = PyList_New(dim1);

    for (i = 0; i < dim1; ++i)
    {
        PyObject* list_i = PyList_New(dim2);
        for (j = 0; j < dim2; ++j)
        {
            PyObject* python_double = Py_BuildValue("d", mat[i][j]);
            PyList_SetItem(list_i, j, python_double);
        }
        PyList_SetItem(returned_list, i, list_i);
    }

    return returned_list;
}

static void matrix_free(matrix *mat, int num_of_rows)
{
    int i;
    for (i=0; i < num_of_rows; i++)
    {
        free((*mat)[i]);
    }
    free((*mat));
}



static PyObject* wam(PyObject *self, PyObject *args)
{
    PyObject *datapoints_obj, *result;
    point *datapoints;
    int N,d;
    matrix wam;
    
    /* parse arguments from python int our variables */
    if (!PyArg_ParseTuple(args, "Oii", &datapoints_obj, &N, &d))
    {
        printf("An Error Has Occurred\n");
        exit(1);  
    }

    /* copy input from python objects into c arrays 
       (read_input returns -1 if function failed) */
    if (read_input(datapoints_obj, &datapoints, N, d) == -1)
    {
        printf("An Error Has Occurred\n");
        exit(1);
    }

    wam = func_wam(datapoints, N, d);
    
    /* convert wam to pyobject */
    result = matrix_to_python(wam, N, N);

    /* free allocated memory */
    matrix_free(&wam, N);
    matrix_free(&datapoints, N);

    return result;
}

static PyObject* ddg(PyObject *self, PyObject *args)
{
    PyObject *wam_obj, *result;
    matrix wam, ddg;
    int N;
    
    /* parse arguments from python int our variables */
    if (!PyArg_ParseTuple(args, "Oi", &wam_obj, &N))
    {
        printf("An Error Has Occurred\n");
        exit(1);  
    }

    /* copy input from python objects into c arrays 
       (read_input returns -1 if function failed) */
    if (read_input(wam_obj, &wam, N, N) == -1)
    {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    
    ddg = func_ddg(wam, N);
    
    /* convert ddg to pyobject */
    result = matrix_to_python(ddg, N, N);

    /* free allocated memory */
    matrix_free(&wam, N);
    matrix_free(&ddg, N);

    return result;
}

static PyObject* lnorm(PyObject *self, PyObject *args)
{
    PyObject *ddg_obj, *wam_obj, *result;
    matrix wam, ddg, lnorm;
    int N;
    
    /* parse arguments from python int our variables */
    if (!PyArg_ParseTuple(args, "OOi", &wam_obj, &ddg_obj, &N))
    {
        printf("An Error Has Occurred\n");
        exit(1);  
    }

    /* copy input from python objects into c arrays 
       (read_input returns -1 if function failed) */
    if ((read_input(wam_obj, &wam, N, N) == -1) || 
        (read_input(ddg_obj, &ddg, N, N) == -1)) 
    {
        if (wam != NULL)
            matrix_free(&wam, N);
        if (ddg != NULL)
            matrix_free(&ddg, N);

        printf("An Error Has Occurred\n");
        exit(1);
    }

    lnorm = func_lnorm(wam, ddg, N);
    
    /* convert lnorm to pyobject */
    result = matrix_to_python(lnorm, N, N);

    /* free allocated memory */
    matrix_free(&wam, N);
    matrix_free(&ddg, N);
    matrix_free(&lnorm, N);
    
    return result;
}

static PyObject* jacobi(PyObject *self, PyObject *args)
{
    PyObject *A_obj, *result;
    matrix A, jacobi;
    int N;
    
    /* parse arguments from python into our variables */
    if (!PyArg_ParseTuple(args, "Oi", &A_obj, &N))
    {
        printf("An Error Has Occurred\n");
        exit(1);  
    }

    /* copy input from python objects into c arrays 
       (read_input returns -1 if function failed) */
    if (read_input(A_obj, &A, N, N) == -1)
    {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    
    jacobi = func_jacobi(A, N);

    result = matrix_to_python(jacobi, N+1, N);

    /* free allocated memory */
    matrix_free(&A, N);
    matrix_free(&jacobi, N+1);

    return result;
}



static int cmpfunc (const void * a, const void * b) 
{
    double *row_a = *((double **)a);
    double *row_b = *((double **)b);
    
    /* return row_a[0] - row_b[0] <= 0 ? -1 : 1 ; */
    if (row_a[0] < row_b[0])
        return -1;
    if (row_b[0] < row_a[0])
        return 1;
    return 0;

}

static void sort_jacobi(matrix jacobi, matrix transpose_jacobi, int N)
{
    int i,j;
    
    /* create transpose_jacobi */
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N+1; j++)
        {
            transpose_jacobi[i][j] = jacobi[j][i];
        }
    }

    qsort(transpose_jacobi, N, sizeof(double *), cmpfunc);
    
    /* converting back to original */
    for (i = 0; i < N+1; i++)
    {
        for (j = 0; j < N; j++)
        {
            jacobi[i][j] = transpose_jacobi[j][i];
        }
    }
}

static int determine_k(matrix jacobi, int N, int k)
{
    /* If k!=0, return as is. Otherwise, perform heuristic */
    if (k != 0)
        return k;

    int argmax, i;
    double max;

    argmax = 0;
    max = jacobi[0][1] - jacobi[0][0];

    for (i = 1; i < floor(N/2); i++)
        if (jacobi[0][i+1] - jacobi[0][i] > max)
        {
            max = jacobi[0][i+1] - jacobi[0][i];
            argmax = i;
        }

    k = argmax;
    /* taking into account that according to instructions,
        the first index should be 1 and not 0           */
    k++;

    return k;
}

static void populate_T(matrix T, matrix jacobi, int N, int k)
{
    int i,j;
    double norm_i;

    /* ASSUMES T is initialized with zeroes */
    for (i = 1; i < N+1; i++)
    {
        /* compute the norm of U_i */
        norm_i = 0;
        for (j = 0; j < k; j++)
            norm_i += pow(jacobi[i][j], 2);
        norm_i = sqrt(norm_i);


        /* populate T_i */
        /* if norm_i == 0, this is a zero vector, 
        so T_i should be a zero vector.       
        otherwise, we should 'normalize' T_i */
        for (j = 0; j < k; j++)
        {
            T[i-1][j] = norm_i == 0 ? 0 : jacobi[i][j] / norm_i;
        }
    }
}

/* perform steps 3-5 of the spectral clustering algorithm */
/* input: jacobi matrix, (N+1) x (N) */
static PyObject* get_input_for_kmeans(PyObject *self, PyObject *args)
{
    PyObject *jacobi_obj, *result;
    matrix jacobi, transpose_jacobi, T;
    int N, k;

    /* parse arguments from python int our variables */
    if (!PyArg_ParseTuple(args, "Oii", &jacobi_obj, &N, &k))
    {
        printf("An Error Has Occurred\n");
        exit(1);  
    }

    /* copy input from python objects into c arrays 
       (read_input returns -1 if function failed) */
    if (read_input(jacobi_obj, &jacobi, N+1, N) == -1)
    {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    
    /* allocate transpose jacobi */
    if (matrix_malloc(&transpose_jacobi, N, N+1) == -1)
    {
        matrix_free(&jacobi, N+1);

        printf("An Error Has Occurred\n");
        exit(1);
    }
    
    
    /* sort jacobi matrix according to the eigenvalues (1st row) */
    sort_jacobi(jacobi, transpose_jacobi, N);

    /* determine new K */
    k = determine_k(jacobi, N, k);

    
    /* https://moodle.tau.ac.il/mod/forum/discuss.php?d=65199 
       (if the heuristic returns k==1, consider this an error) */
    if (k == 1)
    {
        matrix_free(&jacobi, N+1);
        matrix_free(&transpose_jacobi, N);

        printf("An Error Has Occurred\n");
        exit(1);
    }

    /* allocate T */
    if (matrix_malloc(&T, N, k) == -1)
    {
        matrix_free(&jacobi, N+1);
        matrix_free(&transpose_jacobi, N);

        printf("An Error Has Occurred\n");
        exit(1);
    }

    /* populate T */
    populate_T(T, jacobi, N, k);


    /* convert T to pyobject */
    result = matrix_to_python(T, N, k);

    /* free allocated memory */
    matrix_free(&jacobi, N+1);
    matrix_free(&transpose_jacobi, N);
    matrix_free(&T, N);

    return result;
}




/*
----------------------------------------------------------------------- 
-----------------------------API CODE----------------------------------
-----------------------------------------------------------------------
*/



static PyMethodDef capiMethods[] = {
    { "fit", (PyCFunction) fit, METH_VARARGS, PyDoc_STR("kmeans algorithm given initial centroids") },
    { "wam", (PyCFunction) wam, METH_VARARGS, PyDoc_STR("compute weighted adjacency matrix") },
    { "ddg", (PyCFunction) ddg, METH_VARARGS, PyDoc_STR("compute diagonal degree matrix") },
    { "lnorm", (PyCFunction) lnorm, METH_VARARGS, PyDoc_STR("compute normed laplacian matrix") },
    { "jacobi", (PyCFunction) jacobi, METH_VARARGS, PyDoc_STR("compute eigenvalues and eigenvectors") },
    { "get_input_for_kmeans", (PyCFunction) get_input_for_kmeans, METH_VARARGS, PyDoc_STR("compute matrix 'T' for kmeans++") },
    {NULL, NULL, 0, NULL}
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "mykmeanssp", NULL, -1, capiMethods
};


PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        printf("An Error Has Occurred\n");
        exit(1);  
    }
    return m;
}