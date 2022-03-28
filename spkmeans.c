#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <stdlib.h>

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
    int i,j;
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

    /* allocate wam */
    if (matrix_malloc(&wam, N, N) == -1)
    {
        matrix_free(&datapoints, N);

        printf("An Error Has Occurred\n");
        exit(1);
    }

    /* populate wam */
    /* ASSUMES wam is initialized with zeroes */
    /* ASSUMES wam should be symmetric */
    for (i = 0; i < N; i++)
    {
        for (j = i; j < N; j++)
        {
            if (i != j)
                wam[i][j] = exp(-0.5 * norm_of_diff(datapoints[i], datapoints[j], d, 1));
            wam[j][i] = wam[i][j];
        }
    }
    
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
    int i,j;
    
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

    /* allocate ddg */
    if (matrix_malloc(&ddg, N, N) == -1)
    {
        matrix_free(&wam, N);

        printf("An Error Has Occurred\n");
        exit(1);
    }

    /* populate ddg */
    /* ASSUMES ddg is initialized with zeroes */
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            ddg[i][i] += wam[i][j];
        }
    }
    
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
    matrix wam, ddg, lnorm, tmp_mat;
    int N;
    int i,j;
    double tmp;
    
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

    /* allocate lnorm */
    if ((matrix_malloc(&lnorm, N, N) == -1) || 
        (matrix_malloc(&tmp_mat, N, N) == -1))
    {
        matrix_free(&wam, N);
        matrix_free(&ddg, N);

        if (lnorm != NULL)
            matrix_free(&lnorm, N);
        if (tmp_mat != NULL)
            matrix_free(&tmp_mat, N);

        printf("An Error Has Occurred\n");
        exit(1);
    }

    /* populate tmp_mat */
    /* ASSUMES tmp_mat is initialized with zeroes */
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
        {
            tmp = ddg[i][i] == 0 ? 0 : (1 / sqrt(ddg[i][i]));
            tmp_mat[i][j] = tmp * wam[i][j];
        }
                

    /* populate lnorm */
    /* ASSUMES lnorm is initialized with zeroes */
    /* ASSUMES lnorm should be symmetric */
    for (i = 0; i < N; i++)
    {
        for (j = i; j < N; j++)
        {
            tmp = ddg[j][j] == 0 ? 0 : (1 / sqrt(ddg[j][j]));
            lnorm[i][j] = tmp_mat[i][j] * tmp;
        
            lnorm[i][j] = i == j ? 1 - lnorm[i][j] : -lnorm[i][j];

            lnorm[j][i] = lnorm[i][j];
        }
    }
    
    /* convert lnorm to pyobject */
    result = matrix_to_python(lnorm, N, N);

    /* free allocated memory */
    matrix_free(&wam, N);
    matrix_free(&ddg, N);
    matrix_free(&lnorm, N);
    matrix_free(&tmp_mat, N);


    return result;
}



static void update_p_info(matrix A, int N, struct rotation_mat_info *p) 
{
    /* initialize p.i, p.j --- assumes N > 1 */
    p->i = 0 ; p->j = 1 ; 
    int x, y;
    double theta, sign_theta, t;

    /* find the index (i,j) of the largest OFF-DIAGONAL absolute value */
    for (x = 0; x < N; x++)
        for (y = 0; y < N; y++)
            if (x != y && fabs(A[x][y]) > fabs(A[p->i][p->j]))
            {
                p->i = x;
                p->j = y;
            }

    /* obtain c,s */
    theta = (A[p->j][p->j] - A[p->i][p->i]) / (2 * A[p->i][p->j]);
    sign_theta = theta >= 0 ? 1 : -1;
    t = sign_theta / (fabs(theta) + sqrt(theta * theta + 1));
    p->c = 1 / (sqrt(t * t + 1));
    p->s = t * p->c;
}

static void update_A_prime(matrix A_prime, matrix A, int N, struct rotation_mat_info p)
{
    int x,y,r;

    /* start out with values of A */
    for (x = 0; x < N; x++)
        for (y = 0; y < N; y++)
            A_prime[x][y] = A[x][y];
        
    /* make necessary modifications to A_prime */
    for (r = 0; r < N; r++)
    {
        if (r != p.i && r != p.j)
        {
            A_prime[r][p.i] = p.c * A[r][p.i] - p.s * A[r][p.j];
            A_prime[p.i][r] = A_prime[r][p.i];

            A_prime[r][p.j] = p.c * A[r][p.j] + p.s * A[r][p.i];
            A_prime[p.j][r] = A_prime[r][p.j];
        }
    }
    
    A_prime[p.i][p.i] = 
        pow(p.c, 2) * A[p.i][p.i] +
        pow(p.s, 2) * A[p.j][p.j] - 
        2 * p.s * p.c * A[p.i][p.j];

    A_prime[p.j][p.j] = 
        pow(p.s, 2) * A[p.i][p.i] +
        pow(p.c, 2) * A[p.j][p.j] + 
        2 * p.s * p.c * A[p.i][p.j];

    A_prime[p.i][p.j] = 0;
    A_prime[p.j][p.i] = 0;
}

static double off(matrix A, int N)
{
    double val = 0;
    int i,j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
         if (j != i)
            val += pow(A[i][j], 2);
    
    return val;
}

static int converges(matrix A, matrix A_prime, int N)
{
    /* https://moodle.tau.ac.il/mod/forum/discuss.php?d=82078 */
    return (off(A, N) - off(A_prime, N) <= pow(10,-5));
}

static void update_V(matrix V, matrix P, matrix tmp, int N, struct rotation_mat_info p_info)
{
    int x,y,k;

    /* populate P according to p_info */
    for (x = 0; x < N; x++)
        for (y = 0; y < N; y++)
        {
            if ((x == p_info.i && y == p_info.i) || (x == p_info.j && y == p_info.j))
                P[x][y] = p_info.c;
            else
                if (x == y)
                    P[x][y] = 1;
                else 
                    if (x == p_info.i && y == p_info.j)
                        P[x][y] = p_info.s;
                    else
                        if (x == p_info.j && y == p_info.i)
                            P[x][y] = -p_info.s;
                        else
                            P[x][y] = 0;
        }
    

    /* tmp = V*P */
    for (x = 0; x < N; x++)
        for (y = 0; y < N; y++)
        {
            tmp[x][y] = 0;
            for (k = 0; k < N; k++)
                tmp[x][y] += V[x][k] * P[k][y];
        }
            
    /* V <-- tmp */
    for (x = 0; x < N; x++)
        for (y = 0; y < N; y++)
            V[x][y] = tmp[x][y];
        
}

static PyObject* jacobi(PyObject *self, PyObject *args)
{
    PyObject *A_obj;
    matrix A, A_prime, P, V, tmp;
    double *eigenvalues;
    int N,i, x, y;
    int is_diagonal;
    struct rotation_mat_info p_info;
    
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
    

    /* ------------------------------------------------- */

    /* make necessary memory allocations for out matrices */
    if (matrix_malloc(&A_prime, N, N) == -1 || 
        matrix_malloc(&P, N, N) == -1 ||
        matrix_malloc(&V, N, N) == -1 ||
        matrix_malloc(&tmp, N, N))
    {
        matrix_free(&A, N);
        if (A_prime != NULL)
            matrix_free(&A_prime, N);
        if (P != NULL)
            matrix_free(&P, N);
        if (V != NULL)
            matrix_free(&V, N);
        if (tmp != NULL)
            matrix_free(&tmp, N);

        printf("An Error Has Occurred\n");
        exit(1);  
    }

    /* allocate a vector of eigenvalues */
    eigenvalues = malloc(sizeof(double) * N);
    if (eigenvalues == NULL)
    {
        matrix_free(&A, N);
        matrix_free(&A_prime, N);
        matrix_free(&P, N);
        matrix_free(&V, N);
        matrix_free(&tmp, N);

        printf("An Error Has Occurred\n");
        exit(1);  
    }


    /* initialize V as identity matrix */
    /* ASSUMES matrix is initialized with zeroes */
    for (x = 0; x < N; x++)
        V[x][x] = 1;


    /* if A is already a diagonal matrix - do not perform the algorithm */
    is_diagonal = 1;
    for (x = 0; x < N; x++)
        for (y = 0; y < N; y++)
            if (x != y && fabs(A[x][y]) > 0)
                is_diagonal = 0;
    

    if (!is_diagonal)
    {
        /* perform algorithm */
        for (i = 0; i < 100; i++)
        {
            update_p_info(A, N, &p_info);
            update_A_prime(A_prime, A, N, p_info);
            update_V(V, P, tmp, N, p_info);
            if (converges(A, A_prime, N) == 1)
                break;

            /* A <-- A_prime */
            for (x = 0; x < N; x++)
                for (y = 0; y < N; y++)
                    A[x][y] = A_prime[x][y];
        }

        /* assign eigenvalues from algorithm's output */
        for (i = 0; i < N; i++)
            eigenvalues[i] = A_prime[i][i];
    }    
    else
    {
        for (i = 0; i < N; i++)
            eigenvalues[i] = A[i][i];
    }
    /* A <-- A_prime
    for (x = 0; x < N; x++)
        for (y = 0; y < N; y++)
            A[x][y] = A_prime[x][y]; */


    
    /* ------------------------------------------------- */

    /* code is based on this stack-overflow thread:
    https://stackoverflow.com/questions/50668981/how-to-return-a-list-of-ints-in-python-c-api-extension-with-pylist/50683462 */
        
    /* create the object to return */
    PyObject* returned_list = PyList_New(N+1);

    /* add the eigenvalues to the returned object */
    PyObject* list_tmp = PyList_New(N);
    for (x = 0; x < N; ++x)
    {
        PyObject* python_double = Py_BuildValue("d", eigenvalues[x]);
        PyList_SetItem(list_tmp, x, python_double);
    }
    PyList_SetItem(returned_list, 0, list_tmp);

    /* add the eigenvectors to the returned object */
    for (x = 1; x < N+1; ++x)
    {
        PyObject* list_i = PyList_New(N);
        for (y = 0; y < N; ++y)
        {
            PyObject* python_double = Py_BuildValue("d", V[x-1][y]);
            PyList_SetItem(list_i, y, python_double);
        }
        PyList_SetItem(returned_list, x, list_i);
    }


    /* free allocated memory */
    matrix_free(&A, N);
    matrix_free(&A_prime, N);
    matrix_free(&P, N);
    matrix_free(&V, N);
    matrix_free(&tmp, N);
    free(eigenvalues);

    return returned_list;
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