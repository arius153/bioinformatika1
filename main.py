from Bio.Seq import Seq
from Bio import SeqIO;

stops = ["TAA", "TAG", "TGA"]
start = "ATG"


def read_fasta_file(filename):
    for seq_record in SeqIO.parse(filename, "fasta"):
        return seq_record


def split_codones(seq):
    codones = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    if len(codones[-1]) < 3:
        codones.pop()
    return codones


def find_start_stop_codones_without_stop_in_between(codoneslist):
    foundcodones = []
    stopaaa = -1
    for idx, codone in enumerate(codoneslist):
        if idx > stopaaa:
            if codone == start:
                codonecouple = codone
                for idx1, codone2 in enumerate(codoneslist[idx + 1:]):
                    codonecouple = codonecouple + codone2
                    if codone2 in stops:
                        foundcodones.append(codonecouple)
                        stopaaa = idx + idx1
                        break
    return foundcodones


print("1.1 Kodonu poros nereversintos")
print(find_start_stop_codones_without_stop_in_between(split_codones(read_fasta_file("bacterial1.fasta").seq)))
# print("1.2 Kodonu poros reversintos")
# print(find_start_stop_codones_without_stop_in_between(
#     split_codones(read_fasta_file("bacterial1.fasta").seq.reverse_complement())))
