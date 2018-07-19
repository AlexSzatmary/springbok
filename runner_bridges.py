#!/usr/bin/env python
import sys
import pickle
import cloudpickle
import os
import stat


original_path = '/Users/szatmaryac/code/springbok'


N_TASKS_PER_NODE = 28
HEADER = '''#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 02:00:00
#SBATCH --ntasks-per-node '''
original_path = '/Users/szatmaryac/code/springbok'


def make_runner_sh(L_job_name):
    runner_main_path = 'runner.sh'
    with open(runner_main_path, 'w') as hout_main:
        hout_main.write('#!/bin/bash\n')
        for i in range(int(len(L_job_name) / N_TASKS_PER_NODE) + 1):
            runner_sub_path = 'runner-' + str(i) + '.sh'
            hout_main.write('sbatch ' + runner_sub_path + '\n')
            if i == int(len(L_job_name) / N_TASKS_PER_NODE):
                jmax = len(L_job_name) % N_TASKS_PER_NODE
            else:
                jmax = N_TASKS_PER_NODE
            with open(runner_sub_path, 'w') as hout_sub:
                hout_sub.write(HEADER + str(jmax) + '\n')
                hout_sub.write('numactl --show\n')
                hout_sub.write('echo CPUs ${CPUs[0]}\n')
                hout_sub.write(
                    "read -r -a CPUs <<< $(numactl --show|grep physcp|"
                    "sed 's/physcpubind://')\n")
                for j in range(jmax):
                    hout_sub.write(
                        'numactl -C ${CPUs[' + str(j) + ']} python3 runner.py '
                        + L_job_name[i * N_TASKS_PER_NODE + j] + '.pkl &\n')
                hout_sub.write('wait\n')
    st = os.stat(runner_main_path)
    os.chmod(runner_main_path, st.st_mode |
             stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def setup(L_job_name):
    make_runner_sh(L_job_name)
#    scratch_path = os.path.dirname(__file__)
#    scratch_path = os.getcwd()
    # os.system('cp -R ' + scratch_path + ' ' +
    #           os.path.join(original_path, 'batch'))


def main(argv=None):
    if argv is None:
        argv = sys.argv

    (fname, ext) = os.path.splitext(argv[1])

    with open(fname + ext, 'rb') as hin:
        model = pickle.load(hin)
    model.run()
    with open(fname + '.run' + ext, 'wb') as hout:
        cloudpickle.dump(model, hout)

    return model


if __name__ == "__main__":
    sys.exit(main())
