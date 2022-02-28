import subprocess
import sys
import ntpath 
import os


doc_string = \
"""
# --------------------------------------------------------- #
# args (in order):                                          #
#   1) directory where spkmeans.py is                       #   
#   2) directory to output results                          #
#   3) (optional): non-default test-case files              #
# --------------------------------------------------------- #
"""


def handle_args():
    # if no args, print doc_string 
    if (len(sys.argv) < 3 or '?' in sys.argv[1] or '-?' in sys.argv[1]):
        print(doc_string)
        exit()
    
    # parse arguments
    try:
        project_dir, output_dir = sys.argv[1:3]
        tester_dir  = os.path.dirname(sys.argv[0])  # uses this script's directory as tester_dir
        cases_file_path = os.path.join(tester_dir, "test_cases.txt")
        if (len(sys.argv) > 3):
            cases_file_path = os.path.join(tester_dir, sys.argv[3])
    except:
        print("exception when trying to parse arguments...")
        print(sys.argv)
        print(doc_string)
        exit()

    return project_dir, output_dir, tester_dir, cases_file_path


def run_cases(project_dir, output_dir, tester_dir, cases):
    cmd_prefix = ['python', os.path.join(project_dir, 'spkmeans.py')]
    if (os.name != 'nt'):
        cmd_prefix[0] = 'python3'

    for i, case in enumerate(cases):
        # organize command string to run
        case_as_list = case.split(" ")
        case_as_list[-1] = os.path.join(tester_dir, case_as_list[-1])
        case_file_basename = ntpath.basename(case_as_list[-1])

        # run case
        print("running case %d..." % (i+1))
        case_result = subprocess.run(cmd_prefix + case_as_list,
                                     stdout=subprocess.PIPE, 
                                     universal_newlines=True).stdout
        
        # write output to file
        print("writing output...")
        output_file_basename = '_'.join(case_as_list[:-1]) + ' ' + case_file_basename
        output_filename = os.path.join(output_dir, output_file_basename)
        f = open(output_filename, 'w+')
        f.write(case_result)
        f.close()


def main():
    project_dir, output_dir, tester_dir, cases_file_path = handle_args()

    with open(cases_file_path, 'r') as f:
        cases = f.read().split('\n')

    run_cases(project_dir, output_dir, tester_dir, cases)




if __name__ == "__main__":
    main()