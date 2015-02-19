import sys

import pysam
import tables
import argparse
import numpy as np

import genome.db


SNP_UNDEF = -1
SNP_TRACK_NAME = "1000genomes/snp_tab"
SNP_INDEX_TRACK_NAME = "1000genomes/snp_index"

MAX_VAL = 255

def create_carray(track, chrom):
    atom = tables.Int32Atom(dflt=0)
    
    zlib_filter = tables.Filters(complevel=1, complib="zlib")
    
    # create CArray for this chromosome
    shape = [chrom.length]
    carray = track.h5f.createCArray(track.h5f.root, chrom.name,
                                    atom, shape, filters=zlib_filter)

    return carray




def get_snp_index_array(snp_table, chrom):
    snp_indices = np.empty(chrom.length, dtype=np.int32)
    snp_indices[:] = SNP_UNDEF

    if snp_table:
        snp_pos = snp_table[:]['pos']
        idx = range(snp_pos.size)
        snp_indices[snp_pos-1] = idx
        

    return snp_indices
        
    
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--assembly", help="genome assembly that reads "
                        "were mapped to (e.g. hg18)", default=None)

    args = parser.parse_args()

    return args
    

        
def main():
    args = parse_args()
    
    # create a database track
    gdb = genome.db.GenomeDB(assembly=args.assembly)

    snp_index_track = gdb.create_track(SNP_INDEX_TRACK_NAME)
    snp_track = gdb.open_track(SNP_TRACK_NAME)

    chromosomes = gdb.get_chromosomes(get_x=False)
    
    for chrom in chromosomes:
        sys.stderr.write("%s\n" % chrom.name)
        
        index_carray = create_carray(snp_index_track, chrom)

        sys.stderr.write("fetching SNPs\n")
        snp_tab = snp_track.h5f.getNode("/%s" % chrom.name)

        sys.stderr.write("writing SNP index\n")
        snp_index_array = get_snp_index_array(snp_tab, chrom)
        
        index_carray[:] = snp_index_array
    
    snp_index_track.close()
    snp_track.close()

    
main()
        
        
    
