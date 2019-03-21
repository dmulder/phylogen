#!/usr/bin/env python3
from subprocess import Popen, PIPE
import sys, argparse, os
from shutil import which
from glob import glob
from tempfile import NamedTemporaryFile
from multiprocessing.pool import ThreadPool
from multiprocessing import cpu_count

def generate_tree(iqtree, fas_file, out_file_prefix, threads):
    p = Popen([iqtree, '-s', fas_file, '-m', 'TEST', '-nt', '%d' % threads, '-pre', out_file_prefix, '-b', '250'], stdout=PIPE, stderr=PIPE)
    return p.wait() == 0

def iqtree_trees(iqtree, fas_dir, output_dir, cores):
    output_subdir = os.path.join(output_dir, os.path.basename(fas_dir))
    if not os.path.exists(output_subdir):
        os.mkdir(output_subdir)

    procs = 1
    threads = 1
    if int(cores) > 1:
        # 1 Process for every 2 cores allocated
        procs = int(int(cores)/2)
        # 2 threads allowed for every iqtree process
        threads = 2

    fas_files = []
    tree_prefixes = []
    g = '*.fas'
    globs = [g, g.title(), g.upper()]
    for g in globs:
        fas_files.extend(glob(os.path.join(fas_dir, g)))
    tp = ThreadPool(procs)
    for f in fas_files:
        out_file_prefix = os.path.join(output_subdir, os.path.splitext(os.path.basename(f))[0])
        tree_prefixes.append(out_file_prefix)
        tp.apply_async(generate_tree, (iqtree, f, out_file_prefix, threads))
    tp.close()
    tp.join()
    return tree_prefixes

def astral_tree(astral, input_file, output_file):
    cmd = [which('java'), '-jar', astral, '-i', input_file]
    if output_file:
        cmd.extend(['-o', output_file])
    p = Popen(cmd)
    return p.wait() == 0

if __name__ == '__main__':
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-s', help='Input directory of alignments in PHYLIP/FASTA/NEXUS/CLUSTAL/MSF format')
    parser.add_argument('--iqtree', help='Path to iqtree executable', default=os.path.join(current_dir, 'iqtree'))
    parser.add_argument('--astral', help='Path to astral jar executable', default=os.path.join(current_dir, 'astral.jar'))
    parser.add_argument('--cores', help='The number of cores to utilize', default=cpu_count())
    parser.add_argument('--outdir', help='A path to store the iqtree outputs', default=current_dir)
    parser.add_argument('--output', help='a filename for storing the output species tree. Defaults to outputting to stdout.', default=None)
    parser.add_argument('fasdir', help='Directory containing files in fas format')

    args = parser.parse_args()
    iqtree = args.iqtree
    astral = args.astral
    fas_dir = args.fasdir
    output_file = args.output
    output_dir = args.outdir

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

    Popen(['%s -version | grep "version"' % iqtree], shell=True).wait()
    Popen(['%s -jar %s 2>&1 | grep "version"' % (which('java'), astral)], shell=True).wait()

    if not os.path.exists(fas_dir):
        print('The specified fas dir "%s" does not exist' % fas_dir)
        exit(3)
    if not os.path.isdir(fas_dir):
        print('The specified fas dir "%s" is not a directory' % fas_dir)
        exit(4)

    all_trees = None
    tree_prefixes = iqtree_trees(iqtree, os.path.abspath(fas_dir), os.path.abspath(output_dir), args.cores)
    tree_file = os.path.abspath(os.path.join(output_dir, '%s*.treefile' % os.path.basename(fas_dir)))
    with open(tree_file, 'w') as t:
        all_trees = t.name
        for tree_prefix in tree_prefixes:
            treefile = glob('%s*.treefile' % tree_prefix)[-1]
            t.write(open(treefile, 'r').read())

    result = astral_tree(astral, all_trees, output_file)

