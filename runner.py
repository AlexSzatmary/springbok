#!/usr/bin/env python
import sys
import pickle
import cloudpickle
import os
import stat


original_path = '/Users/szatmaryac/code/springbok'


def make_runner_sh(L_job_name):
    runner_sh_path = 'runner.sh'
    with open(runner_sh_path, 'w') as hout:
        hout.write('#!/bin/bash\n')
        hout.write("parallel -j3 --linebuffer ./runner.py {}.pkl ::: " +
                   ' '.join(L_job_name) + '\n')
    st = os.stat(runner_sh_path)
    os.chmod(runner_sh_path, st.st_mode | 
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
