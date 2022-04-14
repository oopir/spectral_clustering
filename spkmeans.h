typedef double* point;

typedef double** matrix;

double norm_of_diff(point v1, point v2, int d, int root);

matrix func_wam(point *datapoints, int N, int d);