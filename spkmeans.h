typedef double* point;

typedef double** matrix;

double norm_of_diff(point v1, point v2, int d, int root);

matrix func_wam(point *datapoints, int N, int d);

matrix func_ddg(matrix wam, int N);

matrix func_lnorm(matrix wam, matrix ddg, int N);

matrix func_jacobi(matrix A, int N);