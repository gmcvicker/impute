import numpy as np
import gzip

def read_geno_vals(geno_file, dtype):
    """Reads haplotypes for the current chromosome into a big matrix"""
    rows = []

    if geno_file.endswith(".gz"):
        f = gzip.open(geno_file)
    else:
        f = open(geno_file)
    
    for line in f:
        words = line.rstrip().split()

        if np.issubdtype(dtype, int):
            vals = np.array([int(x) for x in words[5:]], dtype=dtype)
        else:
            vals = np.array([float(x) for x in words[5:]], dtype=dtype)
        rows.append(vals)

    f.close()
    
    return np.array(rows)

    


# hap_file = "/KG/gmcvicker/IMPUTE/output/CEU/hg19/chr22.hg19.impute2_haps.gz"
# geno_file = "/KG/gmcvicker/IMPUTE/output/CEU/hg19/chr22.hg19.impute2.gz"

hap_file = "/KG/gmcvicker/IMPUTE/output/YRI/hg18/chr22.hg18.impute2_haps.gz"
geno_file = "/KG/gmcvicker/IMPUTE/output/YRI/hg18/chr22.hg18.impute2.gz"

hap_matrix = read_geno_vals(hap_file, dtype=np.uint8)

geno_matrix = read_geno_vals(geno_file, dtype=np.float32)


# ind_idx = 89
ind_idx = 0

hap1 = hap_matrix[:,2*ind_idx]
hap2 = hap_matrix[:,2*ind_idx+1]

homo1_prob = geno_matrix[:,3*ind_idx]
het_prob = geno_matrix[:,3*ind_idx+1]
homo2_prob = geno_matrix[:,3*ind_idx+2]


for i in range(10000):
    #if het_prob[i] > 0.95:
    #if hap1[i] != hap2[i]:
    print "%d %d %.2f %.2f %.2f" %  \
      (hap1[i], hap2[i], homo1_prob[i], het_prob[i], homo2_prob[i])


np.sum(hap1 != hap2)

np.sum(het_prob > 0.95)
np.sum(hap1[(het_prob > 0.95)] != hap2[(het_prob > 0.95)])


np.sum(hap_matrix[:,0] == hap_matrix[:,2*89]).astype(np.float32) / hap_matrix.shape[0]

np.sum(hap_matrix[:,2] == hap_matrix[:,2*89+1]).astype(np.float32) / hap_matrix.shape[0]

np.sum(hap_matrix[:,0] == hap_matrix[:,1]).astype(np.float32) / hap_matrix.shape[0]
