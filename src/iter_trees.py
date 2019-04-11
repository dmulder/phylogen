#!/usr/bin/env python3
import argparse
from glob import glob
import os
from tempfile import NamedTemporaryFile
from shutil import which
from subprocess import Popen, PIPE
import re
# minimum python version 3.5

def astral_tree(astral, input_file, output_file):
    cmd = [which('java'), '-jar', astral, '-i', input_file]
    if output_file:
        cmd.extend(['-o', output_file])
    p = Popen(cmd, stderr=PIPE, stdout=PIPE)
    return p.wait() == 0

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

if __name__ == '__main__':
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('fas_dir', help='Location of fas files used for generating treefiles')
    parser.add_argument('tree_dir', help='Location of tree files generated using fas files in <fas_dir>')

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

    best_rank = 0
    output_dir = './'
    output_file = 'tmp.treefile'
    for i in range(len(tree_files)):
        try:
            os.remove(output_file)
        except:
            pass
        with NamedTemporaryFile('w', dir=output_dir) as t:
            for treefile in tree_files[i:]:
                t.write(open(treefile, 'r').read())

            astral_tree('./astral.jar', t.name, output_file)
            rank = rank_treefile(output_file)
            if rank > best_rank:
                print('Found a better tree using exons from %d to %d' % (i, len(tree_files)))
                best_rank = rank
            else:
                print('Range %d to %d did not improve tree rank' % (i, len(tree_files)))
    print('Best rank was %f' % best_rank)
