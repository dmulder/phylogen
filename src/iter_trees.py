#!/usr/bin/env python3
import argparse
from glob import glob
import os
from tempfile import NamedTemporaryFile
from shutil import which
from subprocess import Popen, PIPE
import re
from multiprocessing import cpu_count
# minimum python version 3.5

current_dir = os.path.dirname(os.path.abspath(__file__))

def astral_tree(astral, input_file, output_file):
    cmd = [which('java'), '-jar', astral, '-i', input_file]
    if output_file:
        cmd.extend(['-o', output_file])
    p = Popen(cmd, stderr=PIPE, stdout=PIPE)
    return p

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
    axon = ''
    with open(fas_file, 'r') as f:
        f.readline()
        axon = f.readline().strip()
    return len(axon)

best_rank = 0
best_output = ''
def validate(p, output_file):
    global best_rank, best_output
    if p.poll() is not None:
        rank = rank_treefile(output_file)
        if rank > best_rank:
            print('Found a better rank of %f in %s' % (rank, output_file))
            best_rank = rank
            best_output = [output_file]
        elif rank == best_rank:
            print('Found the same rank of %f in %s' % (rank, output_file))
            best_output.append(best_output)
        else:
            print('Found a worse rank of %f in %s' % (rank, output_file))
        return True
    return False

if __name__ == '__main__':
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('fas_dir', help='Location of fas files used for generating treefiles')
    parser.add_argument('tree_dir', help='Location of tree files generated using fas files in <fas_dir>')
    parser.add_argument('--cores', help='The number of cores to utilize', default=cpu_count())
    parser.add_argument('--astral', help='Path to astral jar executable', default=os.path.join(current_dir, 'astral.jar'))

    args = parser.parse_args()

    fas_files = []
    g = '*.fas'
    globs = [g, g.title(), g.upper()]
    for g in globs:
        fas_files.extend(glob(os.path.join(args.fas_dir, '**', g), recursive=True))
    if len(fas_files) == 0:
        print('FAS files not found in %s' % args.fas_dir)
        exit(1)

    tree_files = glob(os.path.join(args.tree_dir, '**', '*.treefile'), recursive=True)
    tree_files.sort(key=lambda k : exon_length(k, args.fas_dir))

    output_dir = current_dir
    procs = int(args.cores)-1
    ps = []
    cleanup_files = []
    for i in range(len(tree_files)):
        while len(ps) == procs:
            for p, output_file in ps:
                if validate(p, output_file):
                    ps.remove(p)
        output_file = 'tmp_%dto%d.treefile' % (i, len(tree_files))
        try:
            os.remove(output_file)
        except:
            pass
        with NamedTemporaryFile('w', dir=output_dir, delete=False) as t:
            for treefile in tree_files[i:]:
                t.write(open(treefile, 'r').read())

            print('Calling astral on trees indexed from %d to %d' % (i, len(tree_files)))
            p = astral_tree(args.astral, t.name, output_file)
            ps.append((p, output_file))
            cleanup_files.append(output_file)
    while len(ps) > 0:
        for p, output_file in ps:
            if validate(p, output_file):
                ps.remove(p)
    for cleanup_file in cleanup_files:
        os.remove(cleanup_file)
    print('Best rank was %f with output in %s' % (best_rank, str(best_output)))
