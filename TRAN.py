from BioinformaticsStronghold import grab_fasta


def transitionTransversion(d):
    ids = [_ for _ in d]
    a = ids[0]
    b = ids[1]
    transitions = 0
    transversions = 0
    for i in range(len(d[a])):
        if d[a][i] != d[b][i]:

            if d[a][i] in ['A', 'G'] and d[b][i] in ['A', 'G']:
                transitions += 1
            elif d[a][i] in ['C', 'T'] and d[b][i] in ['C', 'T']:
                transitions += 1
            else:
                transversions += 1
    print(transitions/transversions)


if __name__ == '__main__':
    d = grab_fasta('rosalind_tran.txt')
    transitionTransversion(d)