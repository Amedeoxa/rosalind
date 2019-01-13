def countingDNA(s):
    print(s.count("A"), s.count("C"), s.count("G"), s.count("T"))


if __name__ == '__main__':
    with open('rosalind_dna.txt') as file:
        s = file.readline()
    countingDNA(s)
