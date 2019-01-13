from BioinformaticsStronghold import grab_fasta, translatingRNAtoProtein


def rnaSplicing(d):  # splc RNA splicing, returns protein
    id = [_ for _ in d]
    ss = id[0]
    for i in id[1:]:
        if d[i] in d[ss]:
            d[ss] = d[ss].replace(d[i], '')
        else:
            print('not a substring')

    print(translatingRNAtoProtein(d[ss]))


if __name__ == '__main__':
    d = grab_fasta('rosalind_splc.txt')
    print(rnaSplicing(d))