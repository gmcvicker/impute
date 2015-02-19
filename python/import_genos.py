
import sys

import gzip

import argparse 
import numpy as np
import tables
import os.path

import genome.db




def write_carray(track, chrom, data_type, data_matrix):
    if data_type == "haplotypes":
        atom = tables.Int8Atom(dflt=-1)
    elif data_type == "geno_probs":
        atom = tables.Float32Atom(dflt=np.nan)
    else:
        raise ValueError("unknown datatype %s" % data_type)

    zlib_filter = tables.Filters(complevel=1, complib="zlib")

    shape = data_matrix.shape
    carray = track.h5f.createCArray(track.h5f.root, chrom.name,
                                    atom, shape, filters=zlib_filter)

    carray[:,:] = data_matrix
    
    


def read_samples(sample_file):
    f = open(sample_file)

    samples = []
    for line in f:
        if line.startswith("sample"):
            # header line
            continue
        
        words = line.split()
        samples.append(words[0])

    f.close()

    return samples




def get_n_samples(filename, cols_per_sample):    
    if filename.endswith(".gz"):
        f = gzip.open(filename)
    else:
        f = open(filename)
    
    first_line = f.readline()
    n_words = len(first_line.split())
    n_samples = (n_words - 5) / cols_per_sample
    f.close()

    return n_samples


def make_snp_lookup(snp_tab):
    snp_lookup = {}
    bad_snps = set([])

    row_num = 0
    for row in snp_tab:
        key =  row['pos']
        if key in snp_lookup:
            # Multiple SNPs at same location
            bad_snps.add(key)
        snp_lookup[key] = row_num
        row_num += 1

    return snp_lookup, bad_snps
                

def read_data_matrix(input_path, snp_tab, data_type, chrom):
    if data_type == "haplotypes":
        cols_per_sample = 2
        dtype = np.int8
        dflt = -1
    elif data_type == "geno_probs":
        cols_per_sample = 3
        dtype = np.float32
        dflt = np.nan
    else:
        raise ValueError("unknown datatype %s" % data_type)
    
    # determine number of SNPs and samples
    n_snps = snp_tab.shape[0]
    
    n_samples = get_n_samples(input_path, cols_per_sample)

    sys.stderr.write("  there are %d individuals and %d SNPs\n" %
                     (n_samples, n_snps))

    sys.stderr.write("making SNP lookup table\n")
    snp_lookup, bad_snps = make_snp_lookup(snp_tab)
    
    # pre-allocate matrix of haplotypes
    data_matrix = np.empty((n_snps, n_samples*cols_per_sample), 
                           dtype=dtype)
    
    data_matrix[:] = dflt
    
    # read haplotypes from file into matrix, one SNP at a time
    if input_path.endswith(".gz"):
        f = gzip.open(input_path, "rb")
    else:
        f = open(input_path)

    sys.stderr.write("reading data\n")

    count = 0
    for line in f:        
        count += 1
        if (count % 10000) == 0:
            sys.stderr.write(".")
        
        words = line.rstrip().split()
        key = int(words[2])

        if key in bad_snps:
            continue

        if key in snp_lookup:
            snp_num = snp_lookup[key]
            if np.issubdtype(dtype, int):
                data_matrix[snp_num, :] = [int(x) for x in words[5:]]
            else:
                data_matrix[snp_num, :] = [float(x) for x in words[5:]]
        else:
            sys.stderr.write("WARNING: SNP %d is in data file "
                             "but not in SNP table\n" % key)

    sys.stderr.write("\n")
    
    f.close()

    return data_matrix

    



def parse_args():
    parser = argparse.ArgumentParser(description="Convert genotype "
                                     "probability or haplotype files "
                                     "output by IMPUTE2 into HDF5 format. "
                                     "Input files are assumeed to be named "
                                     "like <input_prefix><chr><input_postfix>")

    parser.add_argument("--assembly", default="hg19", 
                        help="assembly")

    parser.add_argument("--snp_track", default="1000genomes/snp_tab",
                        help="name of track containing SNP information")
    
    parser.add_argument("data_type", choices=("haplotypes", "geno_probs"),
                       help="type of data to import")
    
    parser.add_argument("dest_track", help="name of track to store "
                        "data in (e.g. impute2/yri_geno_probs)")

    parser.add_argument("input_prefix", help="path to input directory "
                        "(e.g. /path/to/dir/)")

    parser.add_argument("--input_postfix", 
                        help="postfix for input files "
                        "(if not specified assumed to be "
                        ".<assembly>.impute2.gz for geno prob files and "
                        ".<assembly>.impute2_haps.gz for haplotype files)",
                        default=None)

    return parser.parse_args()


    

def main():

    args = parse_args()

    gdb = genome.db.GenomeDB(assembly=args.assembly)
        
    snp_track = gdb.open_track(args.snp_track)    
    data_track = gdb.create_track(args.dest_track)

    if args.input_postfix is None:
        if args.data_type == "haplotypes":
            postfix = ".%s.impute2_haps.gz" % args.assembly
        else:
            postfix = ".%s.impute2.gz" % args.assembly
    else:
        postfix = args.input_postfix
        
                
    for chrom in gdb.get_all_chromosomes():
        if not snp_track.has_chromosome(chrom):
            # we don't have SNPs for this chromosome
            continue

        input_filename = args.input_prefix + chrom.name + postfix

        if not os.path.exists(input_filename):
            sys.stderr.write("WARNING skipping %s, input file "
                             "does not exist:\n  %s\n" % (chrom.name, 
                                                        input_filename))
            continue
        
        snp_tab = snp_track.h5f.getNode("/%s" % chrom.name)
        
        data_matrix = read_data_matrix(input_filename, snp_tab, 
                                       args.data_type, chrom)

        sys.stderr.write("writing %s\n" % args.data_type)
        write_carray(data_track, chrom, args.data_type, data_matrix)
    
    data_track.h5f.flush()
    data_track.close()
    snp_track.close()
        

main()
        
