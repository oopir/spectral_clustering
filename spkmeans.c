#define PY_SSIZE_T_CLEAN
#include <Python.h>

typedef double* point;

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


/*  functions that are not exposed */

/* populate an array of points based on a python-object */
/* input: py_obj is a 1d-list, representing a 2d-list with dimensions dim1 x dim2 */
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

double norm_of_diff(point v1, point v2, int d, int root)
{
    int i;
    double result = 0;

    for (i = 0; i < d; i++)
    {
        result += pow((v1[i] - v2[i]), 2);
    }
    
    return root ? sqrt(result) : result;
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


/* free mu & prev_mu:
    - if allocation failed, the failed allocation was on index "**_fail_index"
    - to free only mu / only prev_mu , use NULL,0 for the other parameters */
void free_mu(point *mu, point *prev_mu, int mu_fail_index, int prev_mu_fail_index)
{
    int i;
    for (i = 0; i < mu_fail_index; i++)
    {
        free(mu[i]);
    }
    for (i = 0; i < prev_mu_fail_index; i++)
    {
        free(prev_mu[i]);
    }
    
    if (mu != NULL)
    {
        free(mu);
    }
    if (prev_mu != NULL)
    {
        free(prev_mu);
    }
}

void free_datapoints(point *datapoints, int fail_index)
{
    int i;
    for (i = 0; i < fail_index; i++)
    {
        free(datapoints[i]);
    }
    free(datapoints);
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




static int do_work(point *datapoints, point *mu, int N, int d, int K, int max_iter, double eps)
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
        printf("An Error Has Occurred");
        free_datapoints(datapoints, info.N);
        return 1;
    }

    for (i = 0; i < argum.K; i++)
    {
        clusters[i].sum = malloc(info.d * sizeof(double));
        if (clusters[i].sum == NULL)
        {            
            printf("An Error Has Occurred");
            free_datapoints(datapoints, info.N);
            free_clusters(clusters, i);
            return 1;
        }
    }

    
    prev_mu = malloc(argum.K * sizeof(point));
    if (prev_mu == NULL)
    {            
        printf("An Error Has Occurred");
        free_datapoints(datapoints, info.N);
        free_clusters(clusters, argum.K);
        free_mu(mu, prev_mu, K, -1);
        return 1;
    }

    for (i = 0; i < argum.K; i++)
    {
        prev_mu[i] = malloc(info.d * sizeof(double));
        if (prev_mu[i] == NULL)
        {            
            printf("An Error Has Occurred");
            free_datapoints(datapoints, info.N);
            free_clusters(clusters, argum.K);
            free_mu(mu, prev_mu, K, i);
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
                mu[i][j] = clusters[i].sum[j] / (double)clusters[i].size;
            }
        }

        curr_iter++;
   }

    /* free memory EXCEPT MU (will be freed in "fit") */
    free_datapoints(datapoints, info.N);
    free_clusters(clusters, argum.K);
    free_mu(NULL, prev_mu, 0, argum.K);

    return 0;
}



/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */



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
        printf("An Error has Occurred");
        return NULL;
    }

    /* copy input from python objects into c arrays 
       (read_input returns -1 if function failed) */
    if (read_input(initial_mu_obj, &initial_mu, K, d) == -1 || 
        read_input(datapoints_obj, &datapoints, N, d) == -1)
    {
        printf("An Error has Occurred");
        exit(1);
    }


    /* perform algorithm */
    do_work(datapoints, initial_mu, N, d, K, max_iter, eps);


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

    free_mu(initial_mu, NULL, K, 0);

    return returned_list;
}



static PyMethodDef capiMethods[] = {
    { "fit", (PyCFunction) fit, METH_VARARGS, PyDoc_STR("kmeans algorithm given initial centroids") },
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
        printf("An Error has Occurred");
        return NULL;
    }
    return m;
}
