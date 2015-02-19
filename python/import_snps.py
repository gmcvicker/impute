
import sys
import gzip
import numpy as np
import argparse

import genome.db
import tables



def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--assembly", help="genome assembly (e.g. hg19)",
                        default=None)

    parser.add_argument("--track", help="name of table to store SNPs in",
                        default="impute2/snps")

    parser.add_argument("input_dir", help="base directory. data files are  "
                        "expected to be located in "
                        "<input_dir>/<chr_name><assembly>.impute2.gz")


    return parser.parse_args()
    
    
    


class SNPDesc(tables.IsDescription):
    name = tables.StringCol(16)
    pos = tables.Int32Col()

    # this will truncate the longest variants...
    allele1 = tables.StringCol(32)
    allele2 = tables.StringCol(32)


    


def load_chromosome(options, chrom, chrom_tab, assembly):
    geno_file = "%s/%s.%s.impute2.gz" % (options.input_dir, chrom.name, assembly)

    f = gzip.open(geno_file, "r")

    row = chrom_tab.row

    count = 0
    
    for line in f:
        count += 1
        if count > 10000:
            count = 0
            sys.stderr.write(".")

        
        words = line.rstrip().split()

        row['name'] = words[1]
        row['pos'] = int(words[2])
        row['allele1'] = words[3]
        row['allele2'] = words[4]
        
        row.append()

    chrom_tab.flush()

    sys.stderr.write("\n")
    
    f.close()

    


def main():    
    options = parse_args()

    gdb = genome.db.GenomeDB(assembly=options.assembly)

    track = gdb.create_track(options.track)
    
    chromosomes = gdb.get_chromosomes(get_x=False)
    
    for chrom in chromosomes:
        sys.stderr.write("%s\n" % chrom.name)
        chrom_tab = track.h5f.createTable("/", chrom.name, SNPDesc, "SNPs")
        load_chromosome(options, chrom, chrom_tab, gdb.assembly)

    track.h5f.flush()
    track.close()


main()
        
