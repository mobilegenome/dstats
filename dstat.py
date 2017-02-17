#!/usr/bin/python


from Bio import AlignIO, SeqIO, Align
from sys import argv
import argparse
from os.path import basename, dirname
from math import factorial


blocksize = 5000000  # 5Mb


def do_abbababa(alignment, anc_sequence):
    i = 0
    chr = alignment[0].id
    n = len(alignment)
    anc_sequence = Align.MultipleSeqAlignment([anc_sequence])
    print chr
    for i in xrange(0, len(alignment[1]), blocksize):
        block = alignment[:, i:i + blocksize]
        anc_block = anc_sequence[:, i:i + blocksize]
        f_ab.write("%s\t%i\t%i" % (chr, i, i + len(block[0]) - 1))
        f_extra.write("%s\t%i\t%i" % (chr, i, i + len(block[0]) - 1))
        c = 0
        for h3 in xrange(n):
            for h2 in xrange(n):
                if h2 == h3:
                    continue
                for h1 in xrange(n):
                    if (h1 == h3) or (h1 >= h2):
                        continue
                    print "combination %i of %i" %(c+1, (factorial(n)/2))
                    c += 1
                    #print h1,h2,h3
                    abba = 0
                    baba = 0
                    bbaa = 0


                    for j in range(len(block[1])):  # iterate over sites in alignments

                        s1 = block[h1, j]
                        s2 = block[h2, j]
                        s3 = block[h3, j]
                        s_anc = anc_block[0, j]
                        # print set([h1, h2, h3, h4])
                        if len({s1, s2, s3, s_anc}) != 2:
                            continue  # site not biallelic
                        badchar = False
                        for site in [s1, s2, s3, s_anc]:  # check for N and ambiguities
                            if site in "NYRKMWSBDHV-":
                                badchar = True
                                break

                            # if site in ["N", "Y", "R", "K", "M", "W", "S", "B", "D", "H", "V", "-"]:
                            #     continue
                        if badchar:
                            continue
                        if (s1 == s2) and s3 == s_anc and s1 != s3 and s2 != s_anc:
                                    bbaa += 1
                        elif s1 != s2 and s3 != s_anc:
                            if s1 == s3 and s2 == s_anc:
                                baba += 1
                            elif s2 == s3 and s1 == s_anc:
                                abba += 1
                        else:
                            continue
                    f_ab.write("\t%i\t%i" % (abba, baba))
                    f_extra.write("\t%i\t%i\t%i" % (abba, baba, bbaa))
        f_ab.write("\n")
        f_extra.write("\n")
    return 1


# load as python dictionaries


def main():
    global f_ab, f_extra
    seqs = {}
    records = []
    fname_list = [basename(fpath) for fpath in options.input_files]

    with open(options.output.replace(".abbababa", ".flist"), "w") as fout:
        fout.write("\n".join(options.input_files))

    for fpath in options.input_files:
        fname = basename(fpath)
        seqs[fname] = SeqIO.index(fpath, "fasta")
        records_per_fasta = seqs.get(fname).keys()
        records.extend([record for record in records_per_fasta])
        print fname

    anc = SeqIO.index(options.anc, "fasta")

    print "\n"

    records = set([str(r) for r in records])
    f_ab = open(options.output, "w")
    f_extra = open(options.extra, "w")
    for record in sorted(records):
        sequences = []
        # min_alignment_length = min([len(seqs.get(seq_key).get(record)) for seq_key in fname_list] +
        #                            [len(anc.get(record))])
        for seq_key in fname_list:
            # print seq_key
            sequences.append(seqs.get(seq_key).get(record))

        min_alignment_length = min([len(sequence) for sequence in sequences] + [len(anc.get(record))])

        per_chr_alignment = Align.MultipleSeqAlignment([sequence[:min_alignment_length] for sequence in sequences])

        do_abbababa(per_chr_alignment, anc.get(record)[:min_alignment_length])


    f_ab.close()
    f_extra.close()

    return 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="This is ABBABABA")
    parser.add_argument('-i', '--input_files', required=True,c help='Input files')
    parser.add_argument('-a', '--anc', required=True, help='ancestral fasta')
    parser.add_argument('-o', '--output', required=True, help='*abbababba outputfile')
    parser.add_argument('-e', '--extra', required=True, help='*abbababbabbaa outputfile')

    options = parser.parse_args()

    main()


