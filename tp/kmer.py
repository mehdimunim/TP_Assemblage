def cut_kmer(seq, kmer_len):
    """ 
    Cut a sequence in kmers
    """
    for i in range(len(seq) - kmer_len + 1):
        yield seq[i: i+kmer_len]


def build_kmer_dict(seq, kmer_len):
    """
    Build a dictionary for kmers in sequence
    """
    dic = {}
    for kmer in cut_kmer(seq, kmer_len):
        if kmer in dic:
            dic[kmer] += 1
        else:
            dic.update({kmer: 1})

    return dic


def main():
    seq = "ATTCGGGGCCA"
    len = 2
    # for i in cut_kmer(seq, 4):
    #    print(i)
    print(build_kmer_dict(seq, len))


main()
