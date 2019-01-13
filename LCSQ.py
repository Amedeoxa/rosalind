from BioinformaticsStronghold import grab_fasta

def findingSharedSplicedMotif(d):
    ids = [_ for _ in d]
    lstss = []
    for id in ids:
        lstss.append(d[id])
    lstss.sort(key=len)
    s = lstss[0]
    t = lstss[1]
    lengths = [[0 for j in range(len(t) + 1)] for i in range(len(s) + 1)]
    for i, x in enumerate(s):
        for j, y in enumerate(t):
            if x == y:
                lengths[i + 1][j + 1] = lengths[i][j] + 1
            else:
                lengths[i + 1][j + 1] = max(lengths[i + 1][j], lengths[i][j + 1])
    motif = ''
    x, y = len(s), len(t)
    while x * y != 0:
        if lengths[x][y] == lengths[x - 1][y]:
            x -= 1
        elif lengths[x][y] == lengths[x][y - 1]:
            y -= 1
        else:
            motif = s[x - 1] + motif
            x -= 1
            y -= 1
    print(motif)


if __name__ == '__main__':
    d = grab_fasta('rosalind_lcsq.txt')
    findingSharedSplicedMotif(d)