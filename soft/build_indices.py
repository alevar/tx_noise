#!/usr/bin/env python

# building indices for the simulation

import os
import sys
import glob
import copy
import random
import argparse
import subprocess

# ./build_indices.py --idx_dir /ccb/salz8-1/avaraby/tx_noise/data/indices --reference /ccb/salz3/avaraby/genomes/human/hg38/hg38_p12_ucsc.no_alts.no_fixs.fa --annotation /ccb/salz8-1/avaraby/tx_noise/data/chess2.2_assembly.gff --hisat /ccb/salz8-1/avaraby/agar_manuscript/soft/AGAR2/external/hisat2/hisat2 --salmon /home/avaraby1/soft/salmon-latest_linux_x86_64/bin/salmon --kallisto /home/avaraby1/soft/kallisto/kallisto --gffread gffread --threads 24

def build_indices(args):
    idx_dir = os.path.abspath(args.idx_dir)+"/"
    if not os.path.exists(idx_dir):
        os.makedirs(idx_dir)
    cmd_real = [args.gffread,
                "-w",idx_dir+"annotation.fasta",
                "-g",args.reference,
                args.annotation]
    subprocess.call(cmd_real)

    # build salmon and hisat2 indices
    f=open(idx_dir+"annotation.exons","wb")
    subprocess.call([args.hisat+"_extract_exons.py",
                     args.annotation],
                     stdout=f)
    f.close()
    f=open(idx_dir+"annotation.ss","wb")
    subprocess.call([args.hisat+"_extract_splice_sites.py",
                     args.annotation],
                     stdout=f)
    f.close()
    subprocess.call([args.hisat+"-build",
                     args.reference,
                     idx_dir+"annotation.hisat",
                     "-p",str(args.threads),
                     "--ss",idx_dir+"annotation.ss",
                     "--exon",idx_dir+"annotation.exons"])

    subprocess.call(["salmon","index",
                     "-t",idx_dir+"annotation.fasta",
                     "-p",str(args.threads),
                     "-i",idx_dir+"annotation.salmon"])

    subprocess.call(["kallisto","index",
                     "-i",idx_dir+"annotation.kallisto",
                     idx_dir+"annotation.fasta"])



def main(args):
    parser = argparse.ArgumentParser(description='''Help Page''')

#===========================================
#===================BUILD===================
#===========================================
    parser.add_argument('--idx_dir',
                        required=True,
                        type=str,
                        help="directory to store indices in")
    parser.add_argument('--reference',
                        required=True,
                        type=str,
                        help="path to the reference genome in the fasta format")
    parser.add_argument('--annotation',
                        required=True,
                        type=str,
                        help="path to the annotation which will be provided for each tool")
    parser.add_argument('--hisat',
                        required=True,
                        type=str,
                        help="path to hisat2 executable. Hisat2 path should terminate in the basename of the executable (not directory)")
    parser.add_argument('--salmon',
                        required=True,
                        type=str,
                        help="path to salmon executable")
    parser.add_argument('--kallisto',
                        required=True,
                        type=str,
                        help="path to kallisto executable")
    parser.add_argument("--gffread",
                        required=True,
                        type=str,
                        help="path to gffread executable")
    parser.add_argument("--threads",
                        required=True,
                        type=int,
                        default=1,
                        help="number of threads permitted to use")

    parser.set_defaults(func=build_indices)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main(sys.argv[1:])
