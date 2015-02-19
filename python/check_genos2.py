import numpy as np
import sys
import genome.db


gdb = genome.db.GenomeDB()


sys.stderr.write("checking YRI genos\n")

geno_track = gdb.open_track("impute2/yri_geno_probs")
hap_track = gdb.open_track("impute2/yri_haplotypes")
snp_track = gdb.open_track("impute2/snps")

chr_name = "chr1"

hap_tab = hap_track.h5f.getNode("/%s" % chr_name)
geno_tab = geno_track.h5f.getNode("/%s" % chr_name)
snp_tab = snp_track.h5f.getNode("/%s" % chr_name)

hap1 = hap_tab[:,0]
hap2 = hap_tab[:,1]
het_probs = geno_tab[:,1]

f = ((hap1 != -1) & (hap2 != -1) & (het_probs > 0.95))


sys.stderr.write("%d out of %d sites on %s with hetp > 0.95 have "
                 "het haplotypes\n" % (np.sum(hap1[f] != hap2[f]),
                                       np.sum(f), chr_name))


f = ((hap1 != -1) & (hap2 != -1) & (het_probs < 0.05))

sys.stderr.write("%d out of %d sites on %s with hetp < 0.05 have "
                 "het haplotypes\n" % (np.sum(hap1[f] != hap2[f]),
                                       np.sum(f), chr_name))


sys.stderr.write("checking CEU genos\n")

geno_track = gdb.open_track("impute2/ceu_geno_probs")
hap_track = gdb.open_track("impute2/ceu_haplotypes")
snp_track = gdb.open_track("impute2/snps")

chr_name = "chr1"

hap_tab = hap_track.h5f.getNode("/%s" % chr_name)
geno_tab = geno_track.h5f.getNode("/%s" % chr_name)
snp_tab = snp_track.h5f.getNode("/%s" % chr_name)

hap1 = hap_tab[:,0]
hap2 = hap_tab[:,1]
het_probs = geno_tab[:,1]

f = ((hap1 != -1) & (hap2 != -1) & (het_probs > 0.95))

sys.stderr.write("%d out of %d sites on %s with hetp > 0.95 have "
                 "het haplotypes\n" % (np.sum(hap1[f] != hap2[f]),
                                       np.sum(f), chr_name))

f = ((hap1 != -1) & (hap2 != -1) & (het_probs < 0.05))
np.sum(hap1[f] != hap2[f])
np.sum(f)

sys.stderr.write("%d out of %d sites on %s with hetp < 0.05 have "
                 "het haplotypes\n" % (np.sum(hap1[f] != hap2[f]),
                                       np.sum(f), chr_name))
