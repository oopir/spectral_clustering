import sys
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


def ex2_read_args():
    # validate input format, exit if invalid
    if (len(sys.argv) != 5 and len(sys.argv) != 6) or \
        not sys.argv[1].isnumeric() or \
        (len(sys.argv) == 6 and not sys.argv[2].isnumeric()) or \
        not is_float(sys.argv[-3]) or \
        not (sys.argv[-1].endswith(".txt") or sys.argv[-1].endswith(".csv")) or \
        not (sys.argv[-2].endswith(".txt") or sys.argv[-2].endswith(".csv")):
        print("Invalid Input!")
        sys.exit(1)

    max_iter = 300
    K = int(sys.argv[1])
    eps = float(sys.argv[-3])
    file_name_1 = sys.argv[-2]
    file_name_2 = sys.argv[-1]
    if (len(sys.argv) == 6):
        max_iter = int(sys.argv[2])
    
    if (K < 2 or eps < 0 or max_iter < 0):
        print("Invalid Input!")
        sys.exit(1)

    return K, max_iter, eps, file_name_1, file_name_2


def ex2_main():
    np.random.seed(0)

    # parse arguments
    K, max_iter, eps, file_name_1, file_name_2 = ex2_read_args()

    # load the input data into a dataframe
    merged_data  = pd.merge(pd.read_csv(file_name_1, header=None), pd.read_csv(file_name_2, header=None), on=0)
    merged_data  = merged_data.sort_values(merged_data.columns[0])
    
    datapoints_w_ind = merged_data.to_numpy()
    datapoints       = merged_data.drop(0, axis=1).to_numpy()

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
    data_str = "\n".join([",".join(["%.4f" % fl for fl in mu_i]) for mu_i in result]) + "\n"
    
    print(data_str)


# ----------------------------------------------------------------- #
# ------------------------ end of ex2 code ------------------------ #
# ----------------------------------------------------------------- #


def read_args():
    # if args are invalid, report it 
    # THIS CODE DOES NOT CHECK (k < N) - this is done in main
    if (len(sys.argv) != 4 or \
        not is_int(sys.argv[1]) or \
        sys.argv[2] not in ["spk", "wam", "ddg", "lnorm", "jacobi"] or \
        not (sys.argv[3].endswith(".txt") or sys.argv[3].endswith(".csv"))):
            print("Invalid Input!")
            sys.exit(1)
        
    K         = int(sys.argv[1])
    goal      = sys.argv[2]
    file_name = sys.argv[3]

    if (K < 0 or K == 1):
        print("Invalid Input!")
        sys.exit(1)

    return K, goal, file_name


def pretty_print_mat(mat):
        data_str = "\n".join([",".join(["%.4f" % fl for fl in mat_i]) for mat_i in mat]) + "\n"
        print(data_str)


def main():
    # read args and load input file
    K, goal, file_name = read_args()
    datapoints = pd.read_csv(file_name, header=None).to_numpy()

    # validate arguments based on data files
    if len(datapoints) == 0 or K >= len(datapoints):
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
        
    else:
        print("Invalid Input!")
        exit(1)
    


if __name__ == "__main__":
    main()