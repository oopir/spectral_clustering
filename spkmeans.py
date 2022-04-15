import sys
import os.path
import pandas as pd
import numpy as np
import ctypes
import mykmeanssp

# ----------------------------------------------------------------- #
# --------------------------- ex2 code ---------------------------- #
# ----------------------------------------------------------------- #


def norm_of_diff(v1, v2):
    if (len(v1) != len(v2)):
        print("An Error Has Occurred")
        sys.exit(1)

    norm = sum([(v1[i] - v2[i]) ** 2 for i in range(len(v1))])
    return norm


def is_float(st):
    try:
        float(st)
        return True
    except ValueError:
        return False

def is_int(st):
    try:
        int(st)
        return True
    except ValueError:
        return False


def kmeans_pp(T):
    np.random.seed(0)

    # parse arguments
    K, max_iter, eps = len(T[0]), 300, 0.0

    T_dataframe = pd.DataFrame(T)
    T_dataframe.insert(loc=0, column="a", value=range(len(T)))
    datapoints_w_ind = T_dataframe.to_numpy()

    datapoints       = np.array(T)

    # validate arguments based on data files
    if len(datapoints) == 0 or K >= len(datapoints):
        print("Invalid Input!")
        exit(1)


    # initialize variables needed for the algorithm
    N = len(datapoints)
    d = len(datapoints[0])
    initial_mu = []
    initial_mu_indices = []

    # ------ perform algorithm #1 - choose the initial centroids ------ #

    # choose the first centroid at random
    first_choice_ind = np.random.choice(N)
    initial_mu.append(datapoints[first_choice_ind])
    initial_mu_indices.append(int(datapoints_w_ind[first_choice_ind][0]))
    i = 1

    while i < K:
        dl_list = []
        pl_list = []

        # calculate an updated probability distribution
        for xl in datapoints:
            # next line assumes that initial_mu contains exactly 'i-1' vectors
            dl = min([norm_of_diff(xl, mu_j) for mu_j in initial_mu])
            dl_list.append(dl)

        pl_list = np.divide(dl_list, sum(dl_list))
        
        i += 1

        # choose the next centroid according to the probability distribution
        choice_ind = np.random.choice(N, p=pl_list)
        initial_mu.append(datapoints[choice_ind])
        initial_mu_indices.append(int(datapoints_w_ind[choice_ind][0]))
        
    
    # output the indices of the chosen initial centroids
    print(",".join([str(i) for i in initial_mu_indices]))
    


    # ------ perform algorithm #2 - call c ext to calculate clusters ------ #
    
    result = mykmeanssp.fit(np.array(initial_mu).flatten().tolist(), datapoints.flatten().tolist(), N, d, K, max_iter, eps)
    
    # format the output of the c function
    pretty_print_mat(result)


# ----------------------------------------------------------------- #
# ------------------------ end of ex2 code ------------------------ #
# ----------------------------------------------------------------- #


def read_args():
    # https://moodle.tau.ac.il/mod/forum/discuss.php?d=65199
    # - we can assume K is always provided
    # - if goal != spk, validation on K is NOT required
    # - if goal == spk  --> (k==0) or (1 < k < N)

    # validate all arguments EXCEPT K
    if (len(sys.argv) != 4 or \
        sys.argv[2] not in ["spk", "wam", "ddg", "lnorm", "jacobi"] or \
        not (sys.argv[3].endswith(".txt") or sys.argv[3].endswith(".csv"))):
            print("Invalid Input!")
            sys.exit(1)    
    
    file_name = sys.argv[-1]
    goal      = sys.argv[-2]
    K         = None

    # if goal is spk, CHECK K is valid & read it
    if (goal == 'spk'):
        if (not is_int(sys.argv[1])):
            print("Invalid Input!")
            sys.exit(1)    

        K = int(sys.argv[1])    

        if (K < 0 or K == 1):
            print("Invalid Input!")
            sys.exit(1)

    return K, goal, file_name


def validate_jacobi_input_file(mat):
    #   validating according to the following instructions:
    #   https://moodle.tau.ac.il/mod/forum/discuss.php?d=72817
    #   (input should be a symmetric matrix) 

    if mat.shape[0] != mat.shape[1]:
        print("Invalid Input!")
        exit(1)

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if mat[i][j] != mat[j][i]:
                print("Invalid Input!")
                exit(1)

    
def pretty_print_mat(mat, is_jacobi=False):
        # first condition was inserted because of this post:
        # https://moodle.tau.ac.il/mod/forum/discuss.php?d=70904
        
        first_line = ",".join(["%.4f" % fl for fl in mat[0]]) + "\n"
        if is_jacobi:
            first_line = first_line.replace("-0.0000", "0.0000")
        
        data_str = first_line + "\n".join([",".join(["%.4f" % fl for fl in mat_i]) for mat_i in mat[1:]]) # + "\n"

        print(data_str)


def main():
    # read args and load input file
    K, goal, file_name = read_args()

    try:
        datapoints = pd.read_csv(file_name, header=None).to_numpy()
    except Exception:
        print("Invalid Input!")
        exit(1)

    # validate arguments based on data files
    if len(datapoints) == 0 or (goal == "spk" and K >= len(datapoints)):
        print("Invalid Input!")
        exit(1)
    
    # initialize variables needed for the algorithm
    N = len(datapoints)
    d = len(datapoints[0])


    # act according to the specified goal
    if goal == "wam":
        wam = mykmeanssp.wam(datapoints.flatten().tolist(), N, d)
        pretty_print_mat(wam)

    elif goal == "ddg":
        wam = mykmeanssp.wam(datapoints.flatten().tolist(), N, d)
        wam_flat = sum(wam, [])
        
        ddg = mykmeanssp.ddg(wam_flat, N)
        pretty_print_mat(ddg)
        
    elif goal == "lnorm":
        wam = mykmeanssp.wam(datapoints.flatten().tolist(), N, d)
        wam_flat = sum(wam, [])
        
        ddg = mykmeanssp.ddg(wam_flat, N)
        ddg_flat = sum(ddg, [])

        lnorm = mykmeanssp.lnorm(wam_flat, ddg_flat, N)
        pretty_print_mat(lnorm)
    
    elif goal == "jacobi":
        validate_jacobi_input_file(datapoints)
        jacobi_output = mykmeanssp.jacobi(datapoints.flatten().tolist(), N)
        pretty_print_mat(jacobi_output, is_jacobi=True)

    elif goal == "spk":
        wam = mykmeanssp.wam(datapoints.flatten().tolist(), N, d)
        wam_flat = sum(wam, [])
        
        ddg = mykmeanssp.ddg(wam_flat, N)
        ddg_flat = sum(ddg, [])

        lnorm = mykmeanssp.lnorm(wam_flat, ddg_flat, N)
        lnorm_flat = sum(lnorm, [])

        jacobi_output = mykmeanssp.jacobi(lnorm_flat, N)
        jacobi_flat = sum(jacobi_output, [])
        
        T = mykmeanssp.get_input_for_kmeans(jacobi_flat, N, K)
        
        kmeans_pp(T)


    else:
        print("Invalid Input!")
        exit(1)
    


if __name__ == "__main__":
    main()