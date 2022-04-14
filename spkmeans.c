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

    return wam;
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
    if (!strcmp(".csv", argv[2] + filename_len - 4) && 
        !strcmp(".txt", argv[2] + filename_len - 4))
    {
        printf("Invalid Input!");
        exit(1);
    }
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

void matrix_print(matrix mat, int dim1, int dim2)
{
    int i,j;
    
    for (i = 0; i < dim1; i++)
    {
        for (j = 0; j < (dim2 - 1); j++)
            printf("%.4f,", mat[i][j]);
        printf("%.4f\n", mat[i][dim2-1]);
    }   
}


int main(int argc, char *argv[])
{
    /* declare variables */
    matrix datapoints, wam;
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
        matrix_print(wam, info.N, info.N);
        matrix_free(&wam, info.N);
    }

    return 0;
}
