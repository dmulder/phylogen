#!/usr/bin/env python3
from subprocess import Popen, PIPE
import sys, argparse, os
from shutil import which
from glob import glob
from tempfile import NamedTemporaryFile
from multiprocessing import cpu_count, Pool
import numpy as np
import re
import jpype
import jpype.imports

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
    tree_files = []
    for tree_prefix in tree_prefixes:
        g = '*.treefile'
        globs = [g, g.title(), g.upper()]
        for g in globs:
            tree_files.extend(glob(os.path.join(tree_prefix, g)))
    return tree_files

def astral_tree(args):
    astral, input_file, output_file = args
    jpype.startJVM(jpype.getDefaultJVMPath(), '-Djava.class.path=%s' % astral)
    jpype.imports.registerDomain('phylonet')
    from phylonet.coalescent import CommandLine
    CommandLine.main(['-i', input_file, '-o', output_file])
    jpype.shutdownJVM()

def rank_treefile(filename):
    data = open(filename, 'r').read().strip()
    supports = [float(d) for d in re.findall(r'(\d+\.*\d*):\d+\.*\d*', data)]
    return sum(supports)/len(supports)

def find_matching_fas_file(treefile, fas_dir):
    base = os.path.splitext(os.path.basename(treefile))[0]
    fas_files = []
    g = '.fas'
    globs = [g, g.title(), g.upper()]
    for g in globs:
        fas_files.extend(glob(os.path.join(fas_dir, '**', '*%s*%s' % (base, g)), recursive=True))
    if len(fas_files) == 1:
        return fas_files[0]
    else:
        print('FAS file not found in %s using glob %s' % (fas_dir, os.path.join(fas_dir, '**', '*%s*%s' % (base, '.fas'))))
        exit(1)

def exon_length(treefile, fas_dir):
    fas_file = find_matching_fas_file(treefile, fas_dir)
    exon = ''
    with open(fas_file, 'r') as f:
        f.readline()
        exon = f.readline().strip()
    return len(exon)

best_rank = 0
best_output = []
def validate(output_file):
    global best_rank, best_output
    rank = rank_treefile(output_file)
    if rank > best_rank:
        print('Found a better rank of %f in %s' % (rank, output_file))
        best_rank = rank
        best_output = [output_file]
    elif rank == best_rank:
        print('Found the same rank of %f in %s' % (rank, output_file))
        best_output.append(output_file)
    else:
        print('Found a worse rank of %f in %s' % (rank, output_file))

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
    parser.add_argument('--skip-trees', help='Skip treefile generation, assume they are already created in the output directory', action='store_true')
    parser.add_argument('--start', help='Index to start from (if restarting job)', default=0)
    parser.add_argument('fasdir', help='Directory containing files in fas format')

    args = parser.parse_args()
    iqtree = args.iqtree
    astral = args.astral
    fas_dir = args.fasdir
    final_output_file = args.output
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

    if not args.skip_trees:
        print('Processing trees on %d nodes with %s cores' % (size, args.cores))
        tree_files = iqtree_trees(iqtree, os.path.abspath(fas_dir), os.path.abspath(output_dir), args.cores, args.cores_per_instance)
    else:
        output_subdir = os.path.join(os.path.abspath(output_dir), os.path.basename(os.path.abspath(fas_dir)))
        tree_files = []
        g = '*.treefile'
        globs = [g, g.title(), g.upper()]
        for g in globs:
            tree_files.extend(glob(os.path.join(output_subdir, '**', g)))
    # Only the main node will process the trees via astral
    if rank == 0:
        tree_files.sort(key=lambda k : exon_length(k, fas_dir))
        input_files = []
        output_files = []
        for i in range(int(args.start), len(tree_files)):
            output_file = os.path.join(current_dir, 'tmp_%dto%d.treefile' % (i, len(tree_files)))
            try:
                os.remove(output_file)
            except FileNotFoundError:
                pass
            with NamedTemporaryFile('w', dir=output_dir, delete=False) as t:
                for treefile in tree_files[i:]:
                    t.write(open(treefile, 'r').read())
                input_files.append(t.name)
                output_files.append(output_file)
        ps = Pool(processes=min(int(args.cores)-1, len(input_files)))
        ps.map(astral_tree, zip([args.astral]*len(input_files), input_files, output_files))
        for output_file in output_files:
            validate(output_file)
        if len(best_output) == 1:
            os.rename(best_output[0], final_output_file)
            print('Best rank was %f with output in %s' % (best_rank, final_output_file))
        else:
            print('Best rank was %f with output in %s' % (best_rank, str(best_output)))
        for cleanup_file in input_files + output_files:
            if cleanup_file in best_output:
                continue
            try:
                os.remove(cleanup_file)
            except FileNotFoundError:
                pass
