#!/usr/bin/env python3
from subprocess import Popen, PIPE
import sys, argparse
from shutil import which

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-s', help='Input directory of alignments in PHYLIP/FASTA/NEXUS/CLUSTAL/MSF format')
    parser.add_argument('--iqtree', help='Path to iqtree executable', default=which('iqtree'))
    parser.add_argument('--astral', help='Path to astral jar executable')

    args = parser.parse_args()

    if not args.iqtree:
        print('Path to iqtree not found.\nEither add iqtree to your path, our specify it\'s location via the --iqtree parameter.')
        exit(1)
    if not args.astral:
        print('Path to astral not found.\nSpecify it\'s location via the --astral parameter.')
        exit(1)

    print(args.iqtree)
    Popen(['iqtree -version | grep "version"'], shell=True).wait()

    print(args.astral)
    Popen(['%s -jar %s 2>&1 | grep "version"' % (which('java'), args.astral)], shell=True).wait()
