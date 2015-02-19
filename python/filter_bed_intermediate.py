import sys
import gzip

from operator import attrgetter

MIN_BLOCK_LEN = 100


class SNP(object):
    def __init__(self, pos, rank, line):
        self.pos = pos
        self.rank = rank
        self.line = line
        


def read_snps(chrom, filename):
    """Read SNPs from file, discard ones on - strand and on other chromosome"""
    
    if filename.endswith(".gz"):
        f = gzip.open(filename)
    else:
        f = open(filename)
    
    snp_list = []

    rank = 0
    n_filtered = 0
    
    for l in f:
        words = l.split()

        new_chrom_name = words[0]

        if new_chrom_name != chrom:
            n_filtered += 1
            continue

        strand = words[5]

        if strand != '+':
            n_filtered += 1
            continue
        
        pos = int(words[2])
        rank += 1

        record = words[6]
        
        snp = SNP(pos, rank, record)

        snp_list.append(snp)

    sys.stderr.write("filtered %d SNPs on wrong chromosome or strand\n"
                     % n_filtered)
        
    f.close()
    
    return snp_list
    
        

    
        
def main():
    if len(sys.argv) != 3:
        sys.stderr.write("usage: %s <chr_name> <lifted_over_haplotypes.bed>\n" % 
                         sys.argv[0])
        exit(2)

    chrom = sys.argv[1]
    filename = sys.argv[2]

    sys.stderr.write("reading SNPs\n")
    snp_list = read_snps(chrom, filename)

    # sort snps
    sys.stderr.write("sorting SNPs\n")
    snp_list = sorted(snp_list, key=attrgetter('pos'))
    
    # keep track of blocks of SNPs with ordering that was consistent 
    # before and after sorting
    cur_snp = snp_list[0]
    cur_block = [cur_snp]
    snp_blocks = [cur_block]

    sys.stderr.write("building contiguous blocks of concordant SNPs\n")
    for snp in snp_list[1:]:
        rank_diff = snp.rank - cur_snp.rank

        if rank_diff == 1:
            cur_block.append(snp)
        else:
            # snp rank has changed
            cur_block = [snp]
            snp_blocks.append(cur_block)

        cur_snp = snp

    # now report blocks of SNPs
    sys.stderr.write("writing concordant SNP blocks\n")
    block_num = 1
    for block in snp_blocks:

        block_len = len(block)
        sys.stderr.write("BLOCK #%d, len=%d\n" % (block_num, block_len))

        if block_len < MIN_BLOCK_LEN:
            sys.stderr.write("discarding block below minimum length\n")
            #for snp in block:
            #    sys.stderr.write("#%d %s" % (snp.rank, snp.line))
        else:
            for snp in block:
                record = snp.line.split(",")

                # replace old position with new postition
                record[2] = str(snp.pos)
                sys.stdout.write(" ".join(record) + "\n")

        block_num += 1



    
main()
