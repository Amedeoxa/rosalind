from BioinformaticsStronghold import reverse

def complementingDNA(rna): ##  rosalind_revc.txt
    s = reverse(rna)
    s = s.replace('A', 'O')
    s = s.replace('T', 'A')
    s = s.replace('C', 'X')
    s = s.replace('G', 'C')
    s = s.replace('X', 'G')
    s = s.replace('O', 'T')
    return s


if __name__ == '__main__':
    with open('rosalind_revc.txt') as f:
        rna = f.readline()
    print(complementingDNA(rna))