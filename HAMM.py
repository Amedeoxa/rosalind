def countingPointMutations(s, t):  # hamming distance
    pm = 0
    for i in range(len(t)):
        if s[i] != t[i]:
            pm += 1
    return pm

if __name__ == '__main__':
    with open('rosalind_hamm.txt') as file:
        s = list(file.readline().strip('\n'))  ## remember the end_of_line
        t = list(file.readline().strip('\n'))
    print(countingPointMutations(s, t))