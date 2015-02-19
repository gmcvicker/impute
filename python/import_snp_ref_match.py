# 
# This program records whether the SNP's reference allele
# matches the reference sequence.
#

import sys

import tables
import argparse
import numpy as np

import genome.db


SNP_TRACK_NAME = "impute2/snps"
SNP_REF_MATCH_TRACK_NAME = "impute2/snp_match_ref"
SEQ_TRACK_NAME = "seq"


def create_carray(track, shape, chrom):
    atom = tables.UInt8Atom(dflt=0)
    
    zlib_filter = tables.Filters(complevel=1, complib="zlib")
    
    # create CArray for this chromosome
    carray = track.h5f.createCArray(track.h5f.root, chrom.name,
                                    atom, shape, filters=zlib_filter)

    return carray




    
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--assembly", help="genome assembly that reads "
                        "were mapped to (e.g. hg18)", default=None)
    
    args = parser.parse_args()

    return args
    



        
def main():
    args = parse_args()
    
    gdb = genome.db.GenomeDB(assembly=args.assembly)

    # create a database track
    ref_match_track = gdb.create_track(SNP_REF_MATCH_TRACK_NAME)

    snp_track = gdb.open_track(SNP_TRACK_NAME)
    seq_track = gdb.open_track(SEQ_TRACK_NAME)

    for chrom in gdb.get_chromosomes(get_x=False):
        sys.stderr.write("%s\n" % chrom.name)

        # fetch SNPs and sequence for this chromosome
        sys.stderr.write("fetching SNPs\n")
        snp_tab = snp_track.h5f.getNode("/%s" % chrom.name)

        sys.stderr.write("fetching sequence\n")
        #seq_tab = seq_track.h5f.getNode("/%s" % chrom.name)
        seq_vals = seq_track.get_nparray(chrom)

        ref_match_carray = create_carray(ref_match_track, snp_tab.shape, chrom)
        ref_match_array = ref_match_carray[:]
        
        i = 0
        for snp in snp_tab:
            if i % 10000 == 0:
                sys.stderr.write(".")

            allele_len = len(snp['allele1'])
            start = snp['pos']-1
            end = start + allele_len
            ref_seq = "".join([chr(x) for x in seq_vals[start:end]])

            if ref_seq == snp['allele1']:
                ref_match_array[i] = 1
            else:
                ref_match_array[i] = 0
            i += 1

        
        sys.stderr.write("\n%d/%d SNPs match reference sequence\n" %
                         (np.sum(ref_match_array), i))

        ref_match_carray[:] = ref_match_array
        sys.stderr.write("\n")

        
    snp_track.close()
    seq_track.close()
    ref_match_track.close()

    
main()
        
        
    
