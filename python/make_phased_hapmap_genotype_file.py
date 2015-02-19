import gzip
import sys
import numpy as np

HAPMAP_DIR="/data/share/HapMap/phasing/r22/2008-12-16"



def read_haplotypes(hap_file):
    """Reads haplotypes for the current chromosome into a big matrix"""
    rows = []

    for line in hap_file:
        vals = np.array([int(x) for x in line.rstrip().split()], 
                        dtype=np.uint8)
        rows.append(vals)

    return np.array(rows)

    
def main():
    if len(sys.argv) != 2:
        sys.stderr.write('usage: %s <chromosome>\n' % sys.argv[0])
        exit(2)
    
    chrom = sys.argv[1]

    legend_path = "%s/genotypes_%s_YRI_r22_nr.b36_fwd_legend.txt.gz" % \
      (HAPMAP_DIR, chrom)
    
    phase_path = "%s/genotypes_%s_YRI_r22_nr.b36_fwd.phase.gz" % \
      (HAPMAP_DIR, chrom)

    legend_f = gzip.open(legend_path, "rb")

    # read haplotypes
    phase_f = gzip.open(phase_path, "rb")    
    sys.stderr.write("reading haplotypes for each individual\n")
    hap_array = read_haplotypes(phase_f)
    n_hap = hap_array.shape[1]
    phase_f.close()
    
    header = legend_f.readline()

    sys.stderr.write("writing haplotypes for each SNP\n")
    i = 0
    for line in legend_f:
        legend = line.rstrip().split()

        # write 5 header columns containing identifier (twice),
        # position of SNP, and 2 alleles
        sys.stdout.write("%s %s %s %s %s " % (legend[0], legend[0],
                                              legend[1], legend[2], 
                                              legend[3]))

        # write all genotypes for this SNP, one for each haplotype
        sys.stdout.write(" ".join(["%d" % x for x in hap_array[0:n_hap, i]])
                         + "\n")

        i += 1
        if (i % 100) == 0:
            sys.stderr.write(".")
        
    
    legend_f.close()

    sys.stderr.write("\ndone")
    


main()
