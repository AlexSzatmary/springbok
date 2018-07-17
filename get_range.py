#!/usr/bin/env python
import glob
import numpy as np
import post_processing
import sys
import os
import pickle


def main(argv=None):
    if argv is None:
        argv = sys.argv

    print(argv[1])
    sys.path.append(argv[1])
    a = []
    names = []
    for g in glob.glob(os.path.join(argv[1], '*.run.pkl')):
        print(g)
        with open(g, 'rb') as hin:
            model = pickle.load(hin)
        model.range = post_processing.get_range_continuous(model, 0.5, 60)
        model.meta.append(model.range)
        a.append(model.meta)
        names.append(model.job_name + '\n')
    np.savetxt(os.path.join(argv[1], 'get_range.txt'), a)
    with open(os.path.join(argv[1], 'get_range_names.txt'), 'w') as hout:
        hout.writelines(names)


if __name__ == "__main__":
    sys.exit(main())
