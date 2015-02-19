"""This program takes files of mapped reads and counts the number of
 reads that match the alternate and reference allele at every SNP
 position, as well as the total number of reads starting at each base.
 
 The provided read count files should already be filtered for indel
 overlap etc.

 SNPs are taken from the /impute2/snps track.
 
 Counts of reads overlapping SNPs are stored in the specified
 ref_track, alt_track, other_track.

 Total counts of reads starting at each position are stored in
 the specified read_counts track
"""

import sys
import os
import gzip
import random

import tables
import argparse
import numpy as np

import genome.db

SNP_UNDEF = -1
SNP_TRACK_NAME = "impute2/snps"
HAPLOTYPE_TRACK_NAME = "impute2/all_haplotypes"
SNP_INDEX_TRACK_NAME = "impute2/snp_index"


SAMPLES_FILE = "/data/share/10_IND/IMPUTE/samples.txt"

MAX_VAL = 255



def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--assembly", help="genome assembly that reads "
                        "were mapped to (e.g. hg18)", default=None)

    parser.add_argument("individual",
                        help="identifier of individual to obtain "
                        "heterozygous SNPs for")
    
    parser.add_argument("ref_as_count_track",
                        help="name of track to store counts of reads "
                        "overlapping heterozygous SNPs that match "
                        "reference allele")

    parser.add_argument("alt_as_count_track", 
                        help="name of track to store counts of reads "
                        "overlapping heterozygous SNPs that match "
                        "alternate allele")

    parser.add_argument("other_as_count_track", 
                        help="name of track to store counts of reads "
                        "overlapping heterozygous SNPS that match "
                        "neither reference nor alternate allele")

    parser.add_argument("read_count_track",
                       help="name of track to store counts of "
                       "all reads in--the start (left end) of the "
                       "reads are used")
    
    parser.add_argument("read_filenames", action="store", nargs="+",
                        help="read file(s) to read data from")

    args = parser.parse_args()

    return args
    


def lookup_individual_index(ind_name):
    """Gets the index of individual that is used 
    to lookup information in the genotype and haplotype tables"""
    f = open(SAMPLES_FILE)

    idx = 0
    for line in f:
        words = line.rstrip().split()

        name = words[0].replace("NA", "")
        if name == ind_name:
            f.close()
            return idx
        
        idx += 1

    raise ValueError("individual %s is not in samples file %s" %
                     ind_name, SAMPLES_FILE)




def create_carray(track, chrom):
    atom = tables.UInt8Atom(dflt=0)
    
    zlib_filter = tables.Filters(complevel=1, complib="zlib")
    
    # create CArray for this chromosome
    shape = [chrom.length]
    carray = track.h5f.createCArray(track.h5f.root, chrom.name,
                                    atom, shape, filters=zlib_filter)

    return carray



def get_carray(track, chrom):
    return track.h5f.getNode("/%s" % chrom)




def is_indel(snp):
    if (len(snp['allele1']) != 1) or (len(snp['allele2'])) != 1:
        return True


def add_read_count(read_words, chrom, ref_array, alt_array, other_array,
                   read_count_array, snp_index_array, snp_tab,
                   hap_tab, ind_idx):

    map_code = read_words[4]
    if map_code != "1":
        raise ValueError("expected only uniquely mapping reads, but got "
                         "map_code %s" % map_code)

    read_seq = read_words[0]
    read_len = len(read_seq)

    start = int(read_words[2])
    end = start + read_len-1
    

    # add to total read count
    if read_count_array[start-1] < MAX_VAL:
        read_count_array[start-1] += 1
    else:
        sys.stderr.write("WARNING read count at position %d "
                         "exceeds max %d\n" % (start-1, MAX_VAL))

    # get offsets within read of any overlapping SNPs
    snp_idx = snp_index_array[start-1:end]
    read_offsets = np.where(snp_idx != SNP_UNDEF)[0]

    # get subset of overlapping SNPs that are heterozygous
    het_read_offsets = []
    for offset in read_offsets:
        haps = hap_tab[snp_idx[offset], (ind_idx*2):(ind_idx*2 + 2)]
        if haps[0] != haps[1]:
            # this is a het
            het_read_offsets.append(offset)

    n_overlap_hets = len(het_read_offsets)

    if n_overlap_hets == 0:
        return

    # choose ONE overlapping het SNP randomly to add counts to
    # because we don't want to count same read multiple times
    offset = het_read_offsets[random.randint(0, n_overlap_hets-1)]
    
    snp = snp_tab[snp_idx[offset]]

    if is_indel(snp):
        # sanity check, reads overlapping indels should have been
        # thrown out already
        raise ValueError("read overlaps indel, but reads overlapping "
                         "indels should already have been filtered")
        return

    snp_pos = snp['pos']

    # compare base in read sequence to SNP alleles
    base = read_seq[offset]

    if base == snp['allele1']:
        # matches reference allele
        if ref_array[snp_pos-1] < MAX_VAL:
            ref_array[snp_pos-1] += 1
        else:
            sys.stderr.write("WARNING ref allele count at position %d "
                             "exceeds max %d\n" % (snp_pos, MAX_VAL))
    elif base == snp['allele2']:
        # matches alternate allele
        if alt_array[snp_pos-1] < MAX_VAL:
            alt_array[snp_pos-1] += 1
        else:
            sys.stderr.write("WARNING alt allele count at position %d "
                             "exceeds max %d\n" % (snp_pos, MAX_VAL))
    else:
        # matches neither
        if other_array[snp_pos-1] < MAX_VAL:
            other_array[snp_pos-1] += 1
        else:
            sys.stderr.write("WARNING other allele count at position %d "
                             "exceeds max %d\n" % (snp_pos, MAX_VAL))
                
        
    


def main():
    args = parse_args()
    
    ind_idx = lookup_individual_index(args.individual)
    
    gdb = genome.db.GenomeDB(assembly=args.assembly)

    # create a database tracks to hold read counts
    ref_count_track = gdb.create_track(args.ref_as_count_track)
    alt_count_track = gdb.create_track(args.alt_as_count_track)
    other_count_track = gdb.create_track(args.other_as_count_track)
    read_count_track = gdb.create_track(args.read_count_track)
    
    output_tracks = [ref_count_track, alt_count_track, other_count_track,
                     read_count_track]

    snp_track = gdb.open_track(SNP_TRACK_NAME)
    snp_index_track = gdb.open_track(SNP_INDEX_TRACK_NAME)
    hap_track = gdb.open_track(HAPLOTYPE_TRACK_NAME)

    chrom_dict = {}

    # initialize every chromosome in output tracks
    for chrom in gdb.get_chromosomes(get_x=False):
        for track in output_tracks:
            create_carray(track, chrom)
        chrom_dict[chrom.name] = chrom

    count = 0
    
    
    for read_filename in args.read_filenames:
        sys.stderr.write("reading from file %s\n" % read_filename)
        f = gzip.open(read_filename)

        chrom = None
        for line in f:
            words = line.rstrip().split()

            chrom_name = words[1]
            if chrom_name not in chrom_dict:
                # skip chromosomes we are not interested in, e.g. chrM
                continue

            if chrom is None or chrom_name != chrom.name:
                # store results from previous chromosome
                if chrom:
                    ref_carray[:] = ref_array
                    alt_carray[:] = alt_array
                    other_carray[:] = other_array
                    read_count_carray[:] = read_count_array
                    sys.stderr.write("\n")
                
                # start a new chromosome
                chrom = chrom_dict[chrom_name]
                sys.stderr.write("%s\n" % chrom.name)

                ref_carray = get_carray(ref_count_track, chrom)
                alt_carray = get_carray(alt_count_track, chrom)
                other_carray = get_carray(other_count_track, chrom)
                read_count_carray = get_carray(read_count_track, chrom)

                # retrieve current counts for this chromosome
                ref_array = ref_carray[:]
                alt_array = alt_carray[:]
                other_array = other_carray[:]
                read_count_array = read_count_carray[:]

                # fetch SNPs and sequence for this chromosome
                sys.stderr.write("fetching SNPs\n")
                snp_tab = snp_track.h5f.getNode("/%s" % chrom.name)
                hap_tab = hap_track.h5f.getNode("/%s" % chrom.name)
                snp_index_array = snp_index_track.get_nparray(chrom)

                sys.stderr.write("counting reads\n")

            add_read_count(words, chrom, ref_array, alt_array, 
                           other_array, read_count_array, snp_index_array,
                           snp_tab, hap_tab, ind_idx)

            count += 1
            if count == 10000:
                sys.stderr.write(".")
                count = 0

        # store results from last chromosome
        if chrom:
            ref_carray[:] = ref_array
            alt_carray[:] = alt_array
            other_carray[:] = other_array
            read_count_carray[:] = read_count_array
            sys.stderr.write("\n")

                
        f.close()
        
    for track in output_tracks:
        track.close()
        
    snp_track.close()
    snp_index_track.close()
    hap_track.close()

    
main()
        
        
    
