#!/usr/bin/env python3
from subprocess import Popen, PIPE
import sys, argparse, os
from shutil import which
from glob import glob
from tempfile import NamedTemporaryFile
import multiprocessing

def generate_tree(iqtree, fas_file, out_file_prefix):
    p = Popen([iqtree, '-s', fas_file, '-m', 'TEST', '-nt', '%d' % multiprocessing.cpu_count(), '-pre', out_file_prefix, '-b', '250'], stdout=PIPE, stderr=PIPE)
    return p.wait() == 0

def generate_trees(iqtree, fas_dir):
    fas_files = []
    tree_prefixes = []
    g = '*.fas'
    globs = [g, g.title(), g.upper()]
    for g in globs:
        fas_files.extend(glob(os.path.join(fas_dir, g)))
    for f in fas_files:
        out_file_prefix = None
        with NamedTemporaryFile('w', delete=False) as w:
            out_file_prefix = w.name
        tree_prefixes.append(out_file_prefix)
        generate_tree(iqtree, f, out_file_prefix)
    return tree_prefixes

if __name__ == '__main__':
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-s', help='Input directory of alignments in PHYLIP/FASTA/NEXUS/CLUSTAL/MSF format')
    parser.add_argument('--iqtree', help='Path to iqtree executable', default=os.path.join(current_dir, 'iqtree'))
    parser.add_argument('--astral', help='Path to astral jar executable', default=os.path.join(current_dir, 'astral.jar'))
    parser.add_argument('fasdir', help='Directory containing files in fas format')

    args = parser.parse_args()
    iqtree = args.iqtree
    astral = args.astral
    fas_dir = args.fasdir

    if not os.path.exists(iqtree):
        iqtree = which('iqtree')
    if not iqtree or not os.path.exists(iqtree):
        print('Path to iqtree not found.\nEither add iqtree to your path, or specify it\'s location via the --iqtree parameter.')
        exit(1)
    if not os.path.exists(astral):
        astral = which('astral.jar')
    if not astral or not os.path.exists(astral):
        print('Path to astral not found.\nEither add astral to your path, or specify it\'s location via the --astral parameter.')
        exit(2)

    Popen(['iqtree -version | grep "version"'], shell=True).wait()
    Popen(['%s -jar %s 2>&1 | grep "version"' % (which('java'), args.astral)], shell=True).wait()

    if not os.path.exists(fas_dir):
        print('The specified fas dir "%s" does not exist' % fas_dir)
        exit(3)
    if not os.path.isdir(fas_dir):
        print('The specified fas dir "%s" is not a directory' % fas_dir)
        exit(4)

    tree_prefixes = generate_trees(iqtree, os.path.abspath(fas_dir))
    for tree_prefix in tree_prefixes:
        treefile = glob('%s*.treefile' % tree_prefix)[-1]

