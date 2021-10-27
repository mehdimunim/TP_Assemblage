def fasta_reader(fasta):
    with open(fasta, "r+") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                yield line.rstrip()


def main():
    for i in fasta_reader(r"C:\Users\Mehdi\GitHub\TP_Assemblage\data\eva71_plus_perfect.fq"):
        print(i)


main()
