from Bio import SeqIO
from Bio.SeqUtils import CodonUsage as CU
import os
import glob
import sys

stops = ["TAA", "TAG", "TGA"]
start = "ATG"
MIN_BP_SIZE = 100


def read_fasta_file(filename):
    for seq_record in SeqIO.parse(filename, "fasta"):
        return seq_record


def split_codones(seq):
    seq = str(seq)
    codones = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    if len(codones[-1]) < 3:
        codones.pop()
    return codones


def find_stop_sart_sequences(codoneslist):
    foundcodones = []
    for idx, codone in enumerate(codoneslist):
        if codone in stops:
            startPos = 0
            for idx2, codone2 in enumerate(codoneslist[idx + 1:]):
                if codone2 == start:
                    startPos = idx + idx2 + 1
                if codone2 in stops:
                    if startPos > 0:
                        foundcodones.append(codoneslist[idx:startPos + 1])
                        break
                    else:
                        break
    return foundcodones


def find_sequences(codonesList):
    foundcodonos = []
    insequence = False
    sequence = []
    for idx, codone in enumerate(codonesList):
        if insequence:
            sequence.append(codone)
            if codone in stops:
                sequence.append(codone)
                foundcodonos.append(sequence)
                sequence = []
                insequence = False
        else:
            if codone == start:
                insequence = True
                sequence.append(start)
    return foundcodonos


def find_frames_shorter_than_100(frames):
    for i, frame in enumerate(frames):
        if len(frame) < MIN_BP_SIZE:
            frames.pop(i)
    return frames


ALL_CODONS = ["TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA",
              "TGG",
              "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA",
              "CGG",
              "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA",
              "AGG",
              "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA",
              "GGG"]


def print_frequencies(codonsList):
    codonListLen = len(codonsList)
    for codon in ALL_CODONS:
        count = codonsList.count(codon)
        print(codon + ": " + str(count / codonListLen))


def returnIndex(filePath):
    myIndex = CU.CodonAdaptationIndex()
    myIndex.generate_index(filePath)
    return myIndex


with open("results.txt", 'a') as f:
    sys.stdout = f
    ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
    filePaths = glob.glob(ROOT_DIR + "\data\*.fasta")
    for entry in filePaths:
        print("Rezultatai failo: " + entry)
        readFile = read_fasta_file(entry)
        reversedReadFile = readFile.reverse_complement()
        sequence1 = readFile.seq[:len(readFile.seq)]
        sequence2 = sequence1[1:]
        sequence3 = sequence1[2:]
        revSequence1 = reversedReadFile.seq
        revSequence2 = revSequence1[1:]
        revSequence3 = revSequence1[2:]
        print("\nPirmas punktas:\n")
        pirmoRemo = find_sequences(split_codones(sequence1))
        print("pirmoRemo: ")
        print(pirmoRemo)
        antroRemo = find_sequences(split_codones(sequence2))
        print("antroRemo (pirmoRemo[1:])")
        print(antroRemo)
        trecioRemo = find_sequences(split_codones(sequence3))
        print("trecioRemo (pirmoRemo[2:])")
        print(trecioRemo)
        pirmoApverstoRemo = find_sequences(split_codones(revSequence1))
        print("pirmoApverstoRemo: ")
        print(pirmoApverstoRemo)
        antroApverstoRemo = find_sequences(split_codones(revSequence2))
        print("antroApverstoRemo (pirmoApverstoRemo[1:]): ")
        print(antroApverstoRemo)
        trecioApverstoRemo = find_sequences(split_codones(revSequence3))
        print("trecioApverstoRemo (pirmoApverstoRemo[2:]): ")
        print(trecioApverstoRemo)
        print("\nAntras punktas:\n")
        print("pirmoRemo: ")
        pirmoRemoStopStart = find_stop_sart_sequences(split_codones(sequence1))
        print("\n\n")
        print(len(pirmoRemoStopStart))
        print("\n\n")
        print(pirmoRemoStopStart)
        print("antroRemo (pirmoRemo[1:])")
        antroRemoStopStart = find_stop_sart_sequences(split_codones(sequence2))
        print(antroRemoStopStart)
        print("trecioRemo (pirmoRemo[2:])")
        trecioRemoStopStart = find_stop_sart_sequences(split_codones(sequence3))
        print(trecioRemoStopStart)
        print("pirmoApverstoRemo: ")
        pirmoApverstoRemoStopStart = find_stop_sart_sequences(split_codones(revSequence1))
        print(pirmoApverstoRemoStopStart)
        print("antroApverstoRemo (pirmoApverstoRemo[1:]): ")
        antroApverstoRemoStopStart = find_stop_sart_sequences(split_codones(revSequence2))
        print(antroApverstoRemoStopStart)
        print("trecioApverstoRemo (pirmoApverstoRemo[2:]): ")
        trecioApverstoRemoStopStart = find_stop_sart_sequences(split_codones(revSequence3))
        print(trecioApverstoRemoStopStart)
        print("\nTrecias punktas:\n")
        print("\nStart-Stop: trumpesni uz 100 fragmentu:\n")
        print("pirmoRemo: ")
        print((find_frames_shorter_than_100(pirmoRemo)))
        print("antroRemo (pirmoRemo[1:])")
        print(find_frames_shorter_than_100(antroRemo))
        print("trecioRemo (pirmoRemo[2:])")
        print(find_frames_shorter_than_100(trecioRemo))
        print("pirmoApverstoRemo: ")
        print(find_frames_shorter_than_100(pirmoApverstoRemo))
        print("antroApverstoRemo (pirmoApverstoRemo[1:]): ")
        print(find_frames_shorter_than_100(antroApverstoRemo))
        print("trecioApverstoRemo (pirmoApverstoRemo[2:]): ")
        print(find_frames_shorter_than_100(trecioApverstoRemo))
        print("\nStop-Start: trumpesni uz 100 fragmentu:\n")
        print("pirmoRemo: ")
        print(find_frames_shorter_than_100(pirmoRemoStopStart))
        print("antroRemo (pirmoRemo[1:])")
        print(find_frames_shorter_than_100(antroRemoStopStart))
        print("trecioRemo (pirmoRemo[2:])")
        print(find_frames_shorter_than_100(trecioRemoStopStart))
        print("pirmoApverstoRemo: ")
        print(find_frames_shorter_than_100(pirmoApverstoRemoStopStart))
        print("antroApverstoRemo (pirmoApverstoRemo[1:]): ")
        print(find_frames_shorter_than_100(antroApverstoRemoStopStart))
        print("trecioApverstoRemo (pirmoApverstoRemo[2:]): ")
        print(find_frames_shorter_than_100(trecioApverstoRemoStopStart))
        print('\nKetvirtas punktas\n')
        print('Neapversto kodonu dazniai:')
        print_frequencies(split_codones(sequence1))
        print("Apversto kodonu dažniai:")
        print_frequencies(split_codones(revSequence1))
