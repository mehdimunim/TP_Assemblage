def fasta_reader(fasta):
    with open(fasta, "r+") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                yield line.rstrip()


def main():
    x = fasta_reader(
        r"C:\Users\Mehdi\GitHub\TP_Assemblage\data\eva71_plus_perfect.fq")

    for i in x:
        print(i)


main()
