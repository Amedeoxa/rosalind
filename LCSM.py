from BioinformaticsStronghold import grab_fasta

def findingSharedMotif(d):
    ids = [_ for _ in d]
    lstss = []
    for i in ids:
        lstss.append(d[i])
    lstss = sorted(lstss, key=len)
    shortest = lstss[0]
    others = lstss[1:]
    n = len(lstss[0])
    motif = ''
    for j in range(n):
        for k in range(j, n):
            found = False
            for seq in others:
                if shortest[j:k] in seq:
                    found = True
                else:
                    found = False
                    break

            if found and len(shortest[j:k])> len(motif):
                motif = shortest[j:k]
                mx = len(motif)
    print(mx)
    print(motif)

if __name__ == '__main__':
    d = grab_fasta('rosalind_lcsm.txt')
    findingSharedMotif(d)