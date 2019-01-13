from BioinformaticsStronghold import grab_fasta, countingPointMutations


def distanceMatrix(d):
    ids = [id for id in d]
    mtx = [[] for _ in ids]
    t = len(d[ids[0]])
    for i in range(len(ids)):
        for j in range(len(ids)):

            mtx[i].append(countingPointMutations(d[ids[i]], d[ids[j]])/t)
    for _ in range(len(mtx)):
        print(*mtx[_])


if __name__ == '__main__':
    d = grab_fasta('rosalind_pdst.txt')
    distanceMatrix(d)