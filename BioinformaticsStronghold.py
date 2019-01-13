from collections import OrderedDict
import itertools


def reverse(s):
    strg = ""
    for i in s:
        strg = i + strg
    return strg


def rotate(arr):
    return zip(*arr[::-1])##rotate arr


def getKeysByValue(dictOfElements, valueToFind): # Get a list of keys from dictionary which has the given value
    listOfKeys = list()
    listOfItems = dictOfElements.items()
    for item in listOfItems:
        if item[1] == valueToFind:
            listOfKeys.append(item[0])
    return listOfKeys


def countingDNA():
    lst = []
    with open('rosalind_dna.txt') as file:
        for line in file:
            lst.append(line.split())
    s = lst[0][0]
    print(s.count("A"), s.count("C"), s.count("G"), s.count("T"))


def transcribingDNAintoRNA(t): ##rosalind_rna
    return t.replace('T', 'U')


def complementingDNA(rna): ##  rosalind_revc.txt
    s = reverse(rna)
    s = s.replace('A', 'O')
    s = s.replace('T', 'A')
    s = s.replace('C', 'X')
    s = s.replace('G', 'C')
    s = s.replace('X', 'G')
    s = s.replace('O', 'T')
    return s


def rabbitsAndReacurrence(n, k=1): ## fibonacci
    memo = {0: 0, 1: 1}
    if not n in memo:
        memo[n] = (rabbitsAndReacurrence(n - 1, k) + rabbitsAndReacurrence(n - 2, k)*k)
    return memo[n]


def grab_fasta(filen): # FASTA format. return dictornay computingGCcontent(rosalind_gc.txt)
    d = OrderedDict()
    d2 = OrderedDict()
    with open(filen) as file:
        lst = file.readlines()
    i = 0
    ID = ''
    while i < len(lst):
        if lst[i][:1] == '>':
            ID = lst[i][1:-1]
            d[ID] = []
            i += 1
        else:
            d[ID].append(lst[i].strip('\n'))
            i += 1
    for i in d:
        d[i] = ''.join(d[i])
        d2[i] = (d[i].count('C') + d[i].count('G'))/(len(d[i]))*100
    return d


def countingPointMutations(s, t): ####hamming distance
    #with open('rosalind_hamm.txt') as file:
    #    s = list(file.readline().strip('\n')) ## remember the end_of_line
    #    t = list(file.readline().strip('\n'))

    pm = 0
    for i in range(len(t)):
        if s[i] != t[i]:
            pm += 1
    return pm


def mendelFirstLaw():
    with open('rosalind_iprb.txt') as file:
        k, m, n = map(int, file.readline().split())
        pop = k + m + n
        p = [
            k * (k - 1),  # AA, AA pairs
            k * m,  # AA, Aa pairs
            k * n,  # AA, aa pairs
            m * k,  # Aa, AA pairs
            m * (m - 1) * 0.75,  # Aa, Aa pairs
            m * n * 0.5,  # Aa, aa pairs
            n * k,  # aa, AA pairs
            n * m * 0.5,  # aa, Aa pairs
            0,  # aa, aa pairs
        ]
        pop = k + m + n
    return sum(p)/pop/(pop - 1)


def translatingRNAtoProtein(dna):  ## rosalind_prot.txt
    prot = ''
    d = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

    ##with open(file) as file:
    ##    for line in file:
    ##        rna += line.strip('\n')

    for i in range(0, len(dna), 3):
        if d[dna[i:i+3]] != '_':
            prot += d[dna[i:i+3]]
        else:
            break
    return prot


def findingTheMotif():
    with open('rosalind_subs.txt') as file:
        s = file.readline().strip('\n')
        t = file.readline().strip('\n')
    pos = []
    i = 0
    print(s, t)
    print(t in s)
    while t in s[i:]:
        pos.append(s.find(t, i)+1)
        print(i)
        print(pos)
        i = pos[-1]+1
    print(*pos)


def consensusAndProfile():
    d = grab_fasta('rosalind_cons.txt')
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
    print(prof)
    print(col)
    print(mtx)
    con = list(zip(*prof[::-1]))
    cons = ''
    for i in range(len(prof[0])):
        cons += nuc[3-con[i].index(max(con[i]))]
    print(cons)
    for i in range(4):
        print(nuc[i]+':', *prof[i])


def fibd(n, k=1):  # rosalind_fibd
    ages = [1] + [0]*(k-1)
    for i in range(n-1):
        ages = [sum(ages[1:])] + ages[:-1]
    return sum(ages)


def rnaSplicing():  # splc RNA splicing, returns protein
    d = grab_fasta('rosalind_splc.txt')

    id = [_ for _ in d]
    ss = id[0]
    print(d)
    print(id)
    print(d[id[0]])
    print(ss)
    for i in id[1:]:
        if d[i] in d[ss]:
            print(d[i])
            d[ss] = d[ss].replace(d[i], '')
            print(d[ss])
        else:
            print('not a substring')

    with open('ss.txt', 'w') as file:
        file.write(d[ss])
    return translatingRNAtoProtein('ss.txt')


def overlapGraph():  # overlap graph for the strings in a directed graph Ok (k = len suffix or prefix)
    d = grab_fasta('rosalind_grph.txt')
    ids = [_ for _ in d]
    k = 3
    graphs = []
    graphid = []
    for i in range(len(ids)):
        for j in range(len(ids)):
            if d[ids[i]][len(d[ids[i]])-k:] == d[ids[j]][:-len(d[ids[j]])+k] and d[ids[i]] != d[ids[j]]:
                graphs.append([d[ids[i]], d[ids[j]]])
                graphid.append([ids[i], ids[j]])
    print(graphs)
    for _ in graphid:
        print(*[*_])


def completingTree():
    lst = []
    edgs = []

    with open('rosalind_tree.txt') as file:
        for line in file:
            lst.append(line.strip('\n'))

    n = int(lst[0])
    for i in lst[1:]:
        edgs.append(i.split(' '))

    print(n - len(edgs)-1)  # mmmm


def distanceMatrix():
    d = grab_fasta('rosalind_pdst.txt')
    print(d)
    ids = [id for id in d]

    for id in ids:
        print(len(d[id]))
    mtx = [[] for _ in ids]
    t = len(d[ids[0]])
    for i in range(len(ids)):
        for j in range(len(ids)):

            mtx[i].append(countingPointMutations(d[ids[i]], d[ids[j]])/t)
    for _ in range(len(mtx)):
        print(*mtx[_])


def enumeratingGeneOrders():
    with open('rosalind_perm.txt') as file:
        num = int(file.readline())

    def fact(n):
        fct = 1
        if n == 0:
            return fct
        else:
            return n * fact(n-1)
    a = fact(num)
    prm= list(itertools.permutations(list(range(1, num+1))))
    print(a)
    for _ in range(a):
        print(*prm[_])


def openReadingFrame():  #
    d = grab_fasta('rosalind_orf.txt')
    ids = [_ for _ in d]
    dna = d[ids[0]]
    inv = complementingDNA(dna)
    s = [dna, inv]
    idx = []  # don't really need indexing
    prt = []
    for seq in s:
        for frame in range(3):
            for i in range(frame, len(seq), 3):
                    if seq[i:i+3] == 'ATG':
                        for j in range(i, len(seq), 3):
                            if seq[j:j+3] in ['TAG', 'TAA', 'TGA']:
                                idx.append((i, 'start', frame))  # only appends start if a stop if found
                                idx.append((j+3, 'stop'))
                                prt.append(translatingRNAtoProtein(seq[i:j+3]))
                                break

    print(dna)
    print(inv)
    print(idx)
    for protein in set(prt):  # set for unique proteins
        print(protein)


def locatingRestrictionSites():
    d = grab_fasta('rosalind_revp.txt')
    ids = [_ for _ in d]
    ss = d[ids[0]]
    zz = reverse(complementingDNA(ss))
    rpldr = []
    idx = []
    for length in range(2, 7):
        for j in range(length, len(ss)-length+1):
            z = reverse(zz[j:j+length])
            print([j - length, j], [j, j + length])
            print(ss[j - length:j], z, zz[j:j + length])
            if ss[j-length:j] == z:

                rpldr.append(ss[j-length:j+length])
                idx.append([j-length+1, length*2])
    print(ss)
    print(zz)
    print(rpldr)
    for _ in idx:
        print(*_)


def findingSpliceMotif():
    d = grab_fasta('rosalind_sseq.txt')
    ids = [_ for _ in d]
    ss = d[ids[0]]
    motif = d[ids[1]]
    sub = []
    j = 0
    for i in motif:
        if sub:
            print(i, ss, ss[j:], j)
            sub.append(ss[j:].find(i)+1+j)
            j = sub[-1]
        else:
            sub.append(ss.find(i)+1)

    print(*sub)

def findingSharedMotif():
    d = grab_fasta('rosalind_lcsm.txt')
    ids = [_ for _ in d]
    lstss = []
    for i in ids:
        lstss.append(d[i])
    lstss = sorted(lstss, key=len)
    shortest = lstss[0]
    others = lstss[1:]
    n = len(lstss[0])
    motif = ''
    for j in range(n):
        for k in range(j, n):
            found = False
            for seq in others:
                if shortest[j:k] in seq:
                    found = True
                else:
                    found = False
                    break

            if found and len(shortest[j:k]) > len(motif):
                motif = shortest[j:k]
                mx = len(motif)
    print(mx)
    print(motif)


def transitionTransversion():
    d = grab_fasta('rosalind_tran.txt')
    ids = [_ for _ in d]
    a = ids[0]
    b = ids[1]
    transitions = 0
    transversions = 0
    for i in range(len(d[a])):
        if d[a][i] != d[b][i]:

            if d[a][i] in ['A', 'G'] and d[b][i] in ['A', 'G']:
                print(d[a][i], d[b][i], 'transistion')
                transitions += 1
            elif d[a][i] in ['C', 'T'] and d[b][i] in ['C', 'T']:
                print(d[a][i], d[b][i],'transistion')
                transitions += 1
            else:
                print(d[a][i], d[b][i], '~~~~~~')
                transversions += 1
    print(transitions, transversions, transitions/transversions,)


def mergingTwoArrays():  # rosalind_mer.txt
    with open('test.txt') as file:
       # n1 = int(file.readline())
       # lst1 = list(map(int, file.readline().split()))
       # n2 = int(file.readline())
       # lst2 = list(map(int, file.readline().split()))
        n = int(file.readline())
        lst = list(map(int, file.readline().split()))
        lst1 = lst[:len(lst)//2]
        lst2 = lst[len(lst)//2:]
        n1 = len(lst1)
        n2 = len(lst2)
    res = [0 for _ in range(n1 + n2)]
    i = 0
    j = 0
    k = 0

    while i < n1 and j < n2:
        if lst1[i] < lst2[j]:
            res[k] = lst1[i]
            k += 1
            i += 1
        else:
            res[k] = lst2[j]
            k += 1
            j += 1

    while i < n1:
        res[k] = lst1[i]
        k += 1
        i += 1

    while j < n2:
        res[k] = lst2[j]
        k += 1
        j += 1
    print(*res)


def insertionSort():
    with open('rosalind_ins.txt') as file:
        lst = file.readlines()
    n = int(lst[0])
    array = []
    for x in lst[1:]:
        array = list(map(int, x.split()))
    swaps = 0
    for i in range(1, len(array)):
        k = i - 1
        curr = array[i]
        while (k >= 0) and (array[k] > curr):
            array[k + 1] = array[k]
            swaps += 1
            k -= 1
        array[k + 1] = curr
    return swaps, array


def binarySearch():
    with open('rosalind_bins.txt') as file:
        n = int(file.readline())
        m = int(file.readline())
        lst = list(map(int, file.readline().split()))
        tgt = list(map(int, file.readline().split()))
    rslt = []

    def find(lst, k):
        lower, upper = 0, len(lst) - 1
        while lower <= upper:
            middle = (lower + upper) // 2
            if lst[middle] == k:
                return middle + 1

            elif lst[middle] < k:
                lower = middle + 1
            else:
                upper = middle - 1
        return -1
    for k in tgt:
        rslt.append(find(lst, k))
    print(*rslt)


def merge_sort(arr):
    if len(arr) > 1:
        mid = len(arr) // 2
        lefthalf = arr[:mid]
        righthalf = arr[mid:]

        merge_sort(lefthalf)
        merge_sort(righthalf)

        i = 0
        j = 0
        k = 0
        while i < len(lefthalf) and j < len(righthalf):
            if lefthalf[i] < righthalf[j]:
                arr[k] = lefthalf[i]
                i = i+1
            else:
                arr[k] = righthalf[j]
                j = j+1
            k = k+1

        while i < len(lefthalf):
            arr[k] = lefthalf[i]
            i = i+1
            k = k+1

        while j < len(righthalf):
            arr[k] = righthalf[j]
            j = j+1
            k = k+1
    return arr


def findingSharedSplicedMotif():
    d = grab_fasta('rosalind_lcsq.txt')
    print(d)
    ids = [_ for _ in d]
    lstss = []
    for id in ids:
        lstss.append(d[id])
    print(lstss)

    lstss.sort(key=len)
    s = lstss[0]
    t = lstss[1]
    lengths = [[0 for j in range(len(t) + 1)] for i in range(len(s) + 1)]
    for i, x in enumerate(s):
        for j, y in enumerate(t):
            if x == y:
                lengths[i + 1][j + 1] = lengths[i][j] + 1
            else:
                lengths[i + 1][j + 1] = max(lengths[i + 1][j], lengths[i][j + 1])

    spliced_motif = ''
    x, y = len(s), len(t)
    while x * y != 0:
        if lengths[x][y] == lengths[x - 1][y]:
            x -= 1
        elif lengths[x][y] == lengths[x][y - 1]:
            y -= 1
        else:
            spliced_motif = s[x - 1] + spliced_motif
            x -= 1
            y -= 1
    print(spliced_motif)
