
import sys

import gzip
import numpy as np

import genome.db

import scipy.stats


CHROM = 'chr22'
IMPUTE_DIR = "/data/share/10_IND/IMPUTE"

OUTPUT_FILE = "/KG/gmcvicker/IMPUTE/hapmap_impute_compare.txt"

SAMPLES_FILE = IMPUTE_DIR + "/samples.txt"

HAPMAP_DIR = "/KG/gmcvicker/IMPUTE/hapmap_haplotypes/phaseII/r22/hg18"


def read_samples():
    f = open(SAMPLES_FILE)

    samples = []
    
    for line in f:
        words = line.split()
        samples.append(words[0])

    f.close()

    return samples
    


def read_hapmap_genos(hapmap_file):
    geno_dict = {}
    allele_dict = {}
    
    f = gzip.open(hapmap_file)

    for line in f:
        words = line.rstrip().split()

        n_samp = (len(words) - 5) / 2

        genos = np.empty(n_samp, dtype=np.int8)

        rs_id = words[0]
        allele1 = words[3]
        allele2 = words[4]
        
        for i in range(n_samp):
            idx = i*2 + 5
            geno_call = int(words[idx]) + int(words[idx+1])
            genos[i] = geno_call

        geno_dict[rs_id] = genos
        allele_dict[rs_id] = (allele1, allele2)
        
    f.close()

    return geno_dict, allele_dict



def read_impute_genos(impute_file, hapmap_genos, hapmap_alleles):
    f = gzip.open(impute_file)

    impute_geno_dict = {}
    
    for line in f:
        words = line.rstrip().split()

        n_samp = (len(words) - 5) / 3

        rs_id = words[1]

        genos = np.empty(n_samp, dtype=np.float32)
        
        # only record if this was also in hapmap
        if rs_id in hapmap_genos:

            allele1 = words[3]
            allele2 = words[4]

            hm_alleles = hapmap_alleles[rs_id]

            flip_alleles = False
            if (allele1 == hm_alleles[0]) and (allele2 == hm_alleles[1]):
                flip_alleles = False
            elif (allele1 == hm_alleles[1]) and (allele2 == hm_alleles[0]):
                flip_alleles = True
            else:
                sys.stderr.write("WARNING: mismatch between hapmap and "
                                 "impute alleles for %s: %s/%s != %s/%s\n" %
                                 (rs_id, allele1, allele2, hm_alleles[0], 
                                  hm_alleles[1]))
            
            for i in range(n_samp):
                idx = i*3 + 5

                # prob0 = float(words[idx])
                prob1 = float(words[idx+1])
                prob2 = float(words[idx+2])                
                geno_call = prob1 + 2.0*prob2

                if flip_alleles:
                    genos[i] = 2.0 - geno_call
                else:
                    genos[i] = geno_call

            impute_geno_dict[rs_id] = genos

    return impute_geno_dict




def compare_genos(samples, hapmap_geno_dict, impute_geno_dict):
    n_snp = len(impute_geno_dict)    
    rs_id_list = impute_geno_dict.keys()

    hapmap_calls = np.empty(n_snp, dtype=np.float32)
    impute_calls = np.empty(n_snp, dtype=np.float32)

    
    for i in range(len(samples)):
        # build genotype call array for first sample
        j = 0
        for rs_id in rs_id_list:
            hapmap_calls[j] = hapmap_geno_dict[rs_id][i]
            j += 1

        best_r = 0.0
        best_match = ""
        
        # build genotype call array for second sample
        for k in range(len(samples)):
            j = 0
            for rs_id in rs_id_list:
                impute_calls[j] = impute_geno_dict[rs_id][k]
                j += 1
        
            (r, p) = scipy.stats.pearsonr(hapmap_calls, impute_calls)
            # sys.stderr.write(" %.3f" % r)
            if r > best_r:
                best_r = r
                best_match = samples[k]

                
        if best_match != samples[i]:
            raise ValueError("best match to %s is %s, R=%.3f\n" %
                             (samples[i], best_match, best_r))

            # sys.stderr.write("\n")

        sys.stderr.write("%s %s %.3f\n" % (samples[i], best_match, best_r))

            

def write_geno_calls(out_f, samples, hapmap_geno_dict, impute_geno_dict):
    n_snp = len(impute_geno_dict)    
    rs_id_list = impute_geno_dict.keys()

    hapmap_calls = np.empty(n_snp, dtype=np.float32)
    impute_calls = np.empty(n_snp, dtype=np.float32)

    # write header
    out_f.write("SNP.ID")
    for samp in samples:
        out_f.write(" %s.HAPMAP %s.IMPUTE" % (samp, samp))
    out_f.write("\n")

    # write calls
    for rs_id in rs_id_list:
        out_f.write(rs_id)

        for i in range(len(samples)):
            out_f.write(" %.2f %.2f" % (hapmap_geno_dict[rs_id][i],
                                        impute_geno_dict[rs_id][i]))

        out_f.write("\n")
    
        
    
    


def main():
    gdb = genome.db.GenomeDB()

    chrom = gdb.get_chromosome('chr22')

    samples = read_samples()
    
    impute_file = "%s/hg18/%s.hg18.impute2.gz" % (IMPUTE_DIR ,chrom.name)

    hapmap_file = "%s/haplotypes_%s.YRI.hg18.txt.gz" % (HAPMAP_DIR, chrom.name)

    sys.stderr.write("reading hapmap genotype calls\n")
    hapmap_geno_dict, hapmap_alleles = read_hapmap_genos(hapmap_file)

    sys.stderr.write("reading imputed genotypes\n")
    impute_geno_dict = read_impute_genos(impute_file, hapmap_geno_dict,
                                         hapmap_alleles)

    sys.stderr.write("hapmap SNPs: %d, shared SNPs: %d\n" % 
                     (len(hapmap_geno_dict), len(impute_geno_dict)))
    

    out_f = open(OUTPUT_FILE, "w")
    sys.stderr.write("outputting genotype calls to file %s\n" % OUTPUT_FILE)
    write_geno_calls(out_f, samples, hapmap_geno_dict, impute_geno_dict)
    out_f.close()


    sys.stderr.write("comparing genotype calls\n")
    compare_genos(samples, hapmap_geno_dict, impute_geno_dict)


main()
    

    
    
    
    
    

    
