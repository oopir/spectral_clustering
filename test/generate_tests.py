from sklearn.datasets import make_blobs
import numpy as np
import sys
import random
import os
import subprocess


doc_string = \
"""
# --------------------------------------------------------- #
# args (in order):                                          #
#   1) action (new_files or new_cases)                      #   
#   2) ???????                                              #
#   3) ???????                                              #
# --------------------------------------------------------- #
"""


def pretty_print_mat(mat):
    return "\n".join([",".join(["%.4f" % fl for fl in mat_i]) for mat_i in mat]) + "\n"        


def generate_input_files(input_files_dir, num):
    for i in range(1,num+1):
        N = random.randint(10, 1000)    
        k = random.randint(2, 10)
        d = random.randint(1, 10)
        
        X, _ = make_blobs(n_samples=N, centers=k, n_features=d, shuffle=True)
        X = np.round(X, 4)

        new_file_path = os.path.join(input_files_dir, "input_{}_n={}_d={}.txt".format(i, N, d))
        with open(new_file_path, "w+") as f:
            f.write(pretty_print_mat(X))
        


def generate_test_cases(input_files_dir, cases_file):
    pass


def main():
    action, input_files_dir = sys.argv[1:-1]
    num = int(sys.argv[-1])


    if action == "new_files":
        generate_input_files(input_files_dir, num)
    elif action == "new_cases":
        pass
    else:
        print("bad arguments: first parameter should be 'new_files' or 'new_cases'")


if __name__ == "__main__":
    main()