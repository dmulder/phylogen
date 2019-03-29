#!/usr/bin/env python3
from subprocess import Popen, PIPE
import sys, argparse, os
from shutil import which
from glob import glob
from tempfile import NamedTemporaryFile
from multiprocessing import cpu_count
import numpy as np

size = 1
rank = 0
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
except ImportError:
    pass

def generate_tree(iqtree, fas_file, out_file_prefix, threads):
    print('Processing %s has begun' % fas_file)
    p = Popen([iqtree, '-quiet', '-s', fas_file, '-m', 'TEST', '-nt', 'AUTO', '-ntmax', '%d' % threads, '-pre', out_file_prefix, '-bb', '1000'])
    return p

def iqtree_trees(iqtree, fas_dir, output_dir, cores, cores_per_instance):
    output_subdir = os.path.join(output_dir, os.path.basename(fas_dir))
    if not os.path.exists(output_subdir):
        os.mkdir(output_subdir)

    if rank == 0:
        fas_files = []
        g = '*.fas'
        globs = [g, g.title(), g.upper()]
        for g in globs:
            fas_files.extend(glob(os.path.join(fas_dir, g)))
        tree_prefixes = [os.path.join(output_subdir, os.path.splitext(os.path.basename(f))[0]) for f in fas_files]
        fas_files_chunks = [list(chunk) for chunk in np.split(np.array(fas_files), size)]
    else:
        tree_prefixes = None
        fas_files_chunks = None

    if size == 1 and rank == 0: # Running on a single node does not require MPI
        fas_files_chunk = fas_files_chunks[0]
    else:
        fas_files_chunk = comm.scatter(fas_files_chunks, root=0)

    procs = int(int(cores)/int(cores_per_instance))
    threads = int(cores_per_instance)

    ps = []
    while len(fas_files_chunk) > 0:
        while len(ps) == procs:
            for p in ps:
                if p.poll() is not None:
                    ps.remove(p)
        f = fas_files_chunk.pop()
        prefix = os.path.splitext(os.path.basename(f))[0]
        out_file_dir = os.path.join(output_subdir, prefix)
        if not os.path.exists(out_file_dir):
            os.mkdir(out_file_dir)
        out_file_prefix = os.path.join(out_file_dir, prefix)
        ps.append(generate_tree(iqtree, f, out_file_prefix, threads))
    for p in ps:
        p.wait()
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
    parser.add_argument('--cores-per-instance', help='The number of cores to utilize per iqtree instance', default=1)
    parser.add_argument('--outdir', help='A path to store the iqtree outputs', default=current_dir)
    parser.add_argument('--output', help='a filename for storing the output species tree.', default=os.path.join(current_dir, 'astral.treefile'))
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

    print('Processing trees on %d nodes with %s cores' % (size, args.cores))
    tree_prefixes = iqtree_trees(iqtree, os.path.abspath(fas_dir), os.path.abspath(output_dir), args.cores, args.cores_per_instance)
    # Only the main node will process the trees via astral
    if rank == 0:
        with NamedTemporaryFile('w', dir=output_dir) as t:
            for tree_prefix in tree_prefixes:
                try:
                    treefile = glob(os.path.join(tree_prefix, '*.treefile'))[-1]
                except IndexError:
                    print('Tree file missing for %s' % tree_prefix)
                    continue
                t.write(open(treefile, 'r').read())

            astral_tree(astral, t.name, output_file)

