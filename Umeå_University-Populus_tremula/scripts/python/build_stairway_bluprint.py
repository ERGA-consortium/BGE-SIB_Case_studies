##!/usr/bin/env python3

import numpy as np
import pandas as pd
import natsort as ns
import allel, zarr, os
from functools import reduce



def arguments():
    import argparse, sys
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputFILE", help="Input .vcf (or .vcf.gz) file.",
    	nargs='?', type=str, required=True)
    parser.add_argument("-r", "--runID", help="ID to name the stairway run",
    	nargs='?', type=str, required=True)
    parser.add_argument("-l", "--chrmLEN", help="Chromosome lenght in bp",
    	nargs='?', type=int, required=True)
    parser.add_argument("-m", "--mutationRATE", help="Assumed mutation rate per site per generation",
    	nargs='?', type=float, required=True)
    parser.add_argument("-g", "--generationTIME", help="Assumed generation time (in years)",
    	nargs='?', type=float, required=True)
    parser.add_argument("-d", "--outDIR", help="Output directory for stairway-plot-v2.",
    	nargs='?', type=str, required=True)
    parser.add_argument("-ds", "--stairwayDIR", help="Installatin directory for stairway-plot-v2.",
    	nargs='?', type=str, required=True)
    parser.add_argument("-o", "--outFILE", help="Output blueprint file for stairway plot v2",
    	nargs='?', type=argparse.FileType('w'), required=True)
    parser.add_argument("-f", "--folded", help="To folded (0) or unfolded (1) the SFS. Default:0",
    	nargs='?', type=int, default=0)
    parser.add_argument("-p", "--trainingPERCENT", help="Percentage of sites for training. Default: 0.67",
    	nargs='?', type=float, default=0.67)
    parser.add_argument("-n", "--numberINPUT", help="Number of input files to be created for each estimation. Default: 200",
    	nargs='?', type=int, default=200)
    parser.add_argument("-b", "--breakPOINTS", help="Number of random break points for each try (separated by white space). Default: 7 15 22 28",
    	nargs='?', type=str, default="7 15 22 28")

    


    ArgP = parser.parse_args()
    return ArgP

def main():

    # get input arguments
    args = arguments()

    stairwat_path = args.stairwayDIR if args.stairwayDIR[-1] == "/" else f'{args.stairwayDIR}/'

    to_current_dir = os.getcwd()
    
    print("loading vcf file")
    callset = allel.read_vcf(args.inputFILE)

    print("getting sample information")
    sample_list = callset['samples']
    n_seq = sample_list.shape[0] * 2

    # print("finding index biallelic SNPs")
    # index_biallelic = np.where(callset["variants/numalt"][:] < 2)[0]

    # print("finding index SNP only")
    # index_snps = np.where(callset["variants/is_snp"][:] == True)[0]

    # print("intersect snp index with biallelic index")
    # index_final = reduce(np.intersect1d,(index_biallelic, index_snps))

    # print("counting allels")
    # allel_count = allel.GenotypeChunkedArray(callset['calldata/GT'][:].take(index_final, axis=0)).count_alleles()

    print("counting allels")
    allel_count = allel.GenotypeChunkedArray(callset['calldata/GT'][:]).count_alleles()
    
    print("making sfs")
    if args.folded == 0:
        sfs_full = allel.sfs_folded(allel_count, n=n_seq)
        fold_state = "true"
    elif args.folded == 1:
        sfs_full = allel.sfs(allel_count[:, 1], n=n_seq)
        fold_state = "false"
    sfs = " ".join([str(i) for i in sfs_full[1:n_seq]]) # make it into n=1 to n-1

    print("writing output file")
    blueprint_text = ["#input setting", 
                      f'popid: {args.runID} # id of the population (no white space)',
                      f'nseq: {n_seq} # number of sequences',
                      f'L: {args.chrmLEN} # total number of observed nucleic sites, including polymorphic and monomorphic',
                      f'whether_folded: {fold_state} # whethr the SFS is folded (true or false)',
                      f'SFS: {sfs} # snp frequency spectrum: number of singleton, number of doubleton, etc. (separated by white space)',
                      f'#smallest_size_of_SFS_bin_used_for_estimation: 1 # default is 1; to ignore singletons, uncomment this line and change this number to 2',
                      f'#largest_size_of_SFS_bin_used_for_estimation: {n_seq-1} # default is n-1; to ignore singletons, uncomment this line and change this number to nseq-2',
                      f'pct_training: {args.trainingPERCENT} # percentage of sites for training',
                      f'nrand: {args.breakPOINTS} # number of random break points for each try (separated by white space)',
                      f'project_dir: {to_current_dir}/{args.outDIR}{args.runID} # project directory',
                      f'stairway_plot_dir: {stairwat_path}stairway_plot_es # directory to the stairway plot files',
                      f'ninput: {args.numberINPUT} # number of input files to be created for each estimation',
                      f'#random_seed: 6\n#output setting',
                      f'mu: {args.mutationRATE} # assumed mutation rate per site per generation',
                      f'year_per_generation: {args.generationTIME} # assumed generation time (in years)',
                      f'#plot setting\nplot_title: {args.runID} # title of the plot'
                      f'xrange: 0.1,10000 # Time (1k year) range; format: xmin,xmax; "0,0" for default',
                      f'yrange: 0,0 # Ne (1k individual) range; format: xmin,xmax; "0,0" for default',
                      f'xspacing: 2 # X axis spacing',
                      f'yspacing: 2 # Y axis spacing',
                      f'fontsize: 12 # Font size']
    output_file = args.outFILE
    output_file.write("\n".join(blueprint_text))



if __name__ == '__main__':
    main()