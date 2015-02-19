import gzip
import sys

import genome.db

import tables
import numpy as np


REF_PANEL_DIR = "/data/share/external_public/1000genomes/release/Phase1/IMPUTE/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono"


MISSING_ALLELE_FREQ = np.nan

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--assembly", help="genome assembly (e.g. hg19)",
                        default=None)

    parser.add_argument("population", help="1000 genomes population (e.g. yri)")

    return parser.parse_args()



def main():
    options = parse_args(assembly=options.assembly)
    
    gdb = genome.db.GenomeDB()
    
    snp_track = gdb.open_track("impute2/snps")

    pop = options.population.lower()
    
    maf_track = gdb.create_track("impute2/%s_minor_allele_freq" % pop)
    aaf_track = gdb.create_track("impute2/%s_alt_allele_freq" % pop)

    chromosomes = gdb.get_chromosomes(get_x=False)

    # read allele frequencies from impute legend files
    sys.stderr.write("reading minor allele freqs\n")

    maf_freq_dict = {}
    aaf_freq_dict = {}
    
    for chrom in chromosomes:
        sys.stderr.write("%s\n" % chrom.name)
        leg_path = "%s/ALL_1000G_phase1integrated_v3_%s_impute.legend.gz" % (
            REF_PANEL_DIR, chrom.name)

        f = gzip.open(leg_path)

        header = f.readline()
        
        for line in f:
            words = line.rstrip().split()
            rs_id = words[0]
            alt_allele_freq = float(words[4])
            minor_allele_freq = float(words[8])
            maf_freq_dict[rs_id] = minor_allele_freq
            aaf_freq_dict[rs_id] = alt_allele_freq

        f.close()

    
    # store allele frequencies in Carrays that are parallel to SNP table
    sys.stderr.write("storing allele freqs\n")
    for chrom in chromosomes:
        sys.stderr.write("  %s\n  " % chrom.name)

        snp_tab = snp_track.h5f.getNode("/%s" % chrom.name)

        maf_carray = create_carray(maf_track, chrom, snp_tab)
        aaf_carray = create_carray(aaf_track, chrom, snp_tab)

        idx = 0
        
        for row in snp_tab:

            if idx % 10000 == 0:
                sys.stderr.write(".")
            
            rs_id = row['name']
            if rs_id in maf_freq_dict:
                maf_carray[idx] = maf_freq_dict[rs_id]
                aaf_carray[idx] = aaf_freq_dict[rs_id]
            else:
                # this snp not present in reference panel
                maf_carray[idx] = MISSING_ALLELE_FREQ
                aaf_carray[idx] = MISSING_ALLELE_FREQ
                                            
            idx += 1
        
        sys.stderr.write("\n")
            
        
        

def create_carray(track, chrom, snp_tab):
    atom = tables.Float32Atom(dflt=0.0)
    
    zlib_filter = tables.Filters(complevel=1, complib="zlib")
    
    # create CArray for this chromosome
    shape = [snp_tab.shape[0]]
    carray = track.h5f.createCArray(track.h5f.root, chrom.name,
                                    atom, shape, filters=zlib_filter)

    return carray





main()
