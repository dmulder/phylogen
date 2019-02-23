#!/usr/bin/env python3
from subprocess import Popen, PIPE
import sys, argparse, os
from shutil import which

if __name__ == '__main__':
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-s', help='Input directory of alignments in PHYLIP/FASTA/NEXUS/CLUSTAL/MSF format')
    parser.add_argument('--iqtree', help='Path to iqtree executable', default=os.path.join(current_dir, 'iqtree'))
    parser.add_argument('--astral', help='Path to astral jar executable', default=os.path.join(current_dir, 'astral.jar'))

    args = parser.parse_args()
    iqtree = args.iqtree
    astral = args.astral

    if not os.path.exists(iqtree):
        iqtree = which('iqtree')
    if not iqtree or not os.path.exists(iqtree):
        print('Path to iqtree not found.\nEither add iqtree to your path, or specify it\'s location via the --iqtree parameter.')
        exit(1)
    if not os.path.exists(astral):
        astral = which('astral.jar')
    if not astral or not os.path.exists(astral):
        print('Path to astral not found.\nEither add astral to your path, or specify it\'s location via the --astral parameter.')
        exit(1)

    Popen(['iqtree -version | grep "version"'], shell=True).wait()
    Popen(['%s -jar %s 2>&1 | grep "version"' % (which('java'), args.astral)], shell=True).wait()
