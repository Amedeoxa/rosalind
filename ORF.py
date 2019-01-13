from BioinformaticsStronghold import grab_fasta, translatingRNAtoProtein, complementingDNA


def openReadingFrame(d):
    ids = [_ for _ in d]
    dna = d[ids[0]]
    inv = complementingDNA(dna)
    s = [dna, inv]
    prt = []
    for seq in s:
        for frame in range(3):
            for i in range(frame, len(seq), 3):
                    if seq[i:i+3] == 'ATG':
                        for j in range(i, len(seq), 3):
                            if seq[j:j+3] in ['TAG', 'TAA', 'TGA']:
                                prt.append(translatingRNAtoProtein(seq[i:j+3]))
                                break

    for protein in set(prt):
        print(protein)


if __name__ == '__main__':
        d = grab_fasta('rosalind_orf.txt')
        openReadingFrame(d)