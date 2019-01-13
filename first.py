def stirngAndLists():
    lst = []
    with open(input('file name')) as file:
        for line in file:
            lst.append(line.split())
    a = lst[0][0]
    b = lst[1]
    print(lst)
    print(a[int(b[0]):int(b[1])+1], a[int(b[2]):int(b[3])+1])\

def conditionsAndLoops():
    lst = []
    with open(input()) as file:
        for line in file:
            lst.append(line.split())
    a, b = map(int, lst[0])
    sum = 0
    count = 0
    for i in range(a, b):
        if i % 2 == 1:
            sum += i
            count += 1
    print(sum, count)

def readingAndWriting():
    lst = []
    with open('rosalind_ini5.txt') as file:
        for line in file:
            lst.append(line.split())
    for i in range(len(lst)):
        if not i % 2:
            print(*lst[i+1])

def introToPythonDictionary():
    lst = []
    with open('rosalind_ini6.txt') as file:
        for line in file:
            lst.append(line.split())
    print(lst[0])
    words = lst[0]
    d = {}
    for i in words:
        if i not in d:
            d[i] = words.count(i)
    print(d)
    for i in d:
        print(i, d[i])
