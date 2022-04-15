#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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

static void matrix_free(matrix *mat, int num_of_rows)
{
    int i;
    for (i=0; i < num_of_rows; i++)
    {
        free((*mat)[i]);
    }
    free((*mat));
}


/*
----------------------------------------------------------------------- 
-----------------------------MAIN CODE---------------------------------
-----------------------------------------------------------------------
*/


matrix func_wam(point *datapoints, int N, int d)
{
    int i,j;
    matrix wam;
    
    /* allocate wam */
    if (matrix_malloc(&wam, N, N) == -1)
    {
        matrix_free(&datapoints, N);

        printf("An Error Has Occurred");
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

    return wam;
}

matrix func_ddg(matrix wam, int N)
{
    matrix ddg;
    int i,j;
    
    /* allocate ddg */
    if (matrix_malloc(&ddg, N, N) == -1)
    {
        matrix_free(&wam, N);

        printf("An Error Has Occurred");
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
    
    return ddg;
}

matrix func_lnorm(matrix wam, matrix ddg, int N)
{
    matrix lnorm, tmp_mat;
    int i,j;
    double tmp;
 
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

        printf("An Error Has Occurred");
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
    
    /* free allocated memory */
    matrix_free(&tmp_mat, N);

    return lnorm;
}



static void update_p_info(matrix A, int N, struct rotation_mat_info *p) 
{
    /* initialize p.i, p.j --- assumes N > 1 */
    int x, y;
    double theta, sign_theta, t;
    p->i = 0 ; p->j = 1 ; 
    

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

matrix func_jacobi(matrix A, int N)
{
    matrix A_prime, P, V, tmp, jacobi;
    double *eigenvalues;
    int i, x, y;
    int is_diagonal;
    struct rotation_mat_info p_info;
    

    /* make necessary memory allocations for out matrices */
    if (matrix_malloc(&A_prime, N, N) == -1 || 
        matrix_malloc(&P, N, N) == -1 ||
        matrix_malloc(&V, N, N) == -1 ||
        matrix_malloc(&tmp, N, N) == -1 ||
        matrix_malloc(&jacobi, N+1, N) == -1)
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
        if (jacobi != NULL)
            matrix_free(&jacobi, N+1);

        printf("An Error Has Occurred");
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
        matrix_free(&jacobi, N+1);

        printf("An Error Has Occurred");
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
    

    /* https://moodle.tau.ac.il/mod/forum/discuss.php?d=95579 */
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


    /* jacobi output will contain V, and not V^T:
       https://moodle.tau.ac.il/mod/forum/discuss.php?d=66683 */
    for (y = 0; y < N; y++)
        jacobi[0][y] = eigenvalues[y];
    for (x = 1; x < N+1; x++)
        for (y = 0; y < N; y++)
            jacobi[x][y] = V[x-1][y];
        

    /* free allocated memory */
    matrix_free(&A_prime, N);
    matrix_free(&P, N);
    matrix_free(&V, N);
    matrix_free(&tmp, N);
    free(eigenvalues);

    return jacobi;
}



/* if the function returns info with N=-1 , it means there was an allocation failure */
struct datapoint_info read_input_file(char* input_filename, point **datapoints)
{
    struct datapoint_info info = {0, 1};
    char ch;
    int i,j;
    double tmp;

    FILE *file = NULL;
    file = fopen(input_filename, "r");
    if (file == NULL)
    {
        printf("Invalid Input!");
        exit(1);   
    }

    /* figure out N and d */
    while((ch=fgetc(file))!=EOF) 
    {
        if (ch=='\n')
        {
            info.N++;
        }
   }
    rewind(file);
    while((ch=fgetc(file))!='\n') 
    {
        if (ch==',')
        {
            info.d++;
        }
   }
    rewind(file);
    
    /* allocate memory to datapoint array */
    *datapoints = malloc(sizeof(point) * info.N);
    if (datapoints == NULL)
    {
        info.N = -1;
        info.d = 0;
        return info;
    }
    
    /* read datapoints */
    for (i = 0; i < info.N; i++)
    {
        (*datapoints)[i] = malloc(sizeof(double) * info.d);
        if ((*datapoints)[i] == NULL)
        {
            info.N = -1;
            info.d = i;
            return info;
        }

        for (j = 0; j < info.d; j++)
        {
            fscanf(file, "%lf", &tmp);
            (*datapoints)[i][j] = tmp;
            fgetc(file);
        }        
    }

    fclose(file);

    return info;
}


void validate_arguments(int argc, char *argv[])
{
    int filename_len;
    
    if (argc != 3  || 
        (strcmp(argv[1],"wam") != 0    &&  strcmp(argv[1],"ddg") != 0 &&  
         strcmp(argv[1],"lnorm") != 0  &&  strcmp(argv[1],"jacobi") != 0))
    {
        printf("Invalid Input!");
        exit(1);
    }
    
    if (strlen(argv[2]) <= 4)
    {
        printf("Invalid Input!");
        exit(1);
    }
    else
        filename_len = strlen(argv[2]);


    /* code to check that string ends with specific extension:
     * https://stackoverflow.com/questions/10347689/how-can-i-
     * check-whether-a-string-ends-with-csv-in-c */
    if (strcmp(".csv", argv[2] + filename_len - 4) && 
        strcmp(".txt", argv[2] + filename_len - 4))
    {
        printf("Invalid Input!");
        exit(1);
    }
}


void validate_jacobi_input_file(matrix datapoints, struct datapoint_info info)
{
    /* validating according to the following instructions:
       https://moodle.tau.ac.il/mod/forum/discuss.php?d=72817
       (input should be a symmetric matrix) */
    
    int i,j;

    if (info.N != info.d)
    {
        printf("Invalid Input!");
        matrix_free(&datapoints, info.N);
        exit(1);
    }
    for (i = 0; i < info.N; i++)
        for (j = 0; j < info.d; j++)
            if (datapoints[i][j] != datapoints[j][i])
            {
                printf("Invalid Input!");
                matrix_free(&datapoints, info.N);
                exit(1);
            }      
}


void matrix_print(matrix mat, int dim1, int dim2, int is_jacobi)
{
    int i,j;
    
    /* first line (if jacobi, make sure -0 is printed as +0) 
       https://moodle.tau.ac.il/mod/forum/discuss.php?d=70904  */
    if (is_jacobi)
    {
        char buffer2[20];

        for (j = 0; j < (dim2 - 1); j++)
        {
            char buffer[20];

            sprintf(buffer, "%.4f", mat[0][j]);
            if (!strcmp(buffer, "-0.0000"))
                printf("0.0000,");
            else
                printf("%.4f,", mat[0][j]);
        }

        sprintf(buffer2, "%.4f", mat[0][dim2-1]);
        if (!strcmp(buffer2, "-0.0000"))
            printf("0.0000\n");
        else
            printf("%.4f\n", mat[0][dim2-1]);
    }
    else
    {    
        for (j = 0; j < (dim2 - 1); j++)
            printf("%.4f,", mat[0][j]);
        printf("%.4f\n", mat[0][dim2-1]);
    }


    /* all other lines */
    for (i = 1; i < dim1; i++)
    {
        for (j = 0; j < (dim2 - 1); j++)
            printf("%.4f,", mat[i][j]);
        printf("%.4f\n", mat[i][dim2-1]);
    }   
}


int main(int argc, char *argv[])
{
    /* declare variables */
    matrix datapoints, wam, ddg, lnorm, jacobi;
    struct datapoint_info info;


    /* read arguments, datapoint info and datapoints */
    validate_arguments(argc, argv);
    info = read_input_file(argv[2], &datapoints);
    /* check if there was an allocation error */
    if (info.N == -1  ||  info.N == 0)
    {
        printf("An Error Has Occurred");
        matrix_free(&datapoints, info.d);
        exit(1);
    }


    if (!strcmp(argv[1],"wam"))
    {
        wam = func_wam(datapoints, info.N, info.d);
        matrix_print(wam, info.N, info.N, 0);
        matrix_free(&wam, info.N);
    }
    if (!strcmp(argv[1],"ddg"))
    {
        wam = func_wam(datapoints, info.N, info.d);
        ddg = func_ddg(wam, info.N);
        
        matrix_print(ddg, info.N, info.N, 0);
        matrix_free(&wam, info.N);
        matrix_free(&ddg, info.N);
    }
    if (!strcmp(argv[1],"lnorm"))
    {
        wam   = func_wam(datapoints, info.N, info.d);
        ddg   = func_ddg(wam, info.N);
        lnorm = func_lnorm(wam, ddg, info.N);

        matrix_print(lnorm, info.N, info.N, 0);
        matrix_free(&wam, info.N);
        matrix_free(&ddg, info.N);
        matrix_free(&lnorm, info.N);
    }
    if (!strcmp(argv[1],"jacobi"))
    {
        validate_jacobi_input_file(datapoints, info);

        jacobi = func_jacobi(datapoints, info.N);
        matrix_print(jacobi, info.N+1, info.N, 1);

        matrix_free(&jacobi, info.N+1);
    }

    matrix_free(&datapoints, info.N);

    return 0;
}
