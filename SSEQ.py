from BioinformaticsStronghold import grab_fasta


def findingSpliceMotif(d):
    ids = [_ for _ in d]
    ss = d[ids[0]]
    motif = d[ids[1]]
    sub = []
    j = 0
    for i in motif:
        if sub:
            sub.append(ss[j:].find(i)+1+j)
            j = sub[-1]
        else:
            sub.append(ss.find(i)+1)

    print(*sub)


if __name__ == '__main__':
    d = grab_fasta('rosalind_sseq.txt')
    findingSpliceMotif(d)