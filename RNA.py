def transcribingDNAintoRNA(t): ##rosalind_rna
    return t.replace('T', 'U')


if __name__ == '__main__':
    with open('rosalind_rna.txt') as f:
        t = f.readline()
    print(transcribingDNAintoRNA(t))