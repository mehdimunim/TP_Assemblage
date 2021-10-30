import os


def save_contigs(contigs_list, output_file):
    with open(output_file, "w+") as out:
        for i, contig in enumerate(contigs_list):
            # header
            out.write(f">contig_{i} len={contig[1]}\n")
            # sequence
            out.write(fill(contig[0])+"\n")


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return "\n".join(text[i: i+width] for i in range(0, len(text), width))


def main():
    file = r"C:\Users\Mehdi\GitHub\TP_Assemblage\tp\test.fastq"
    list = [("AATTC", 10), ("AA"*100, 10), ("TTCCC"*5, 100)]
    #save_contigs(list, file)
    list = [1, 1, 2]
    print(list[None:None])


main()
