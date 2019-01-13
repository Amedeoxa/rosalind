def findingTheMotif(s, t):
    pos = []
    i = 0
    while t in s[i:]:
        pos.append(s.find(t, i)+1)
        i = pos[-1]+1
    print(*pos)


if __name__ == '__main__':
    with open('rosalind_subs.txt') as file:
        s = file.readline().strip('\n')
        t = file.readline().strip('\n')
    findingTheMotif(s, t)
