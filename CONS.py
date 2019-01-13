from BioinformaticsStronghold import grab_fasta


def consensusAndProfile(d):  # fasta dictionary
    mtx = []
    nuc = ['A', 'C', 'G', 'T']
    for i in d:
        mtx.append(list(d[i]))
    col = list(zip(*mtx[::-1]))
    prof = []
    for i in range(4):
        prof.append([])
        for j in col:
            prof[i].append(j.count(nuc[i]))
    con = list(zip(*prof[::-1]))
    cons = ''
    for i in range(len(prof[0])):
        cons += nuc[3-con[i].index(max(con[i]))]
    for i in range(4):
        print(nuc[i]+':', *prof[i])


if __name__ == '__main__':
    d = grab_fasta('rosalind_cons.txt')
    consensusAndProfile(d)