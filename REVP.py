from BioinformaticsStronghold import grab_fasta, reverse, complementingDNA


def locatingRestrictionSites(d):
    ids = [_ for _ in d]
    ss = d[ids[0]]
    zz = reverse(complementingDNA(ss))
    rpldr = []
    idx = []
    for length in range(2, 7):
        for j in range(length, len(ss)-length+1):
            z = reverse(zz[j:j+length])
            if ss[j-length:j] == z:
                rpldr.append(ss[j-length:j+length])
                idx.append([j-length+1, length*2])

    for _ in idx:
        print(*_)


if __name__ == '__main__':
    d = grab_fasta('rosalind_revp.txt')
    locatingRestrictionSites(d)