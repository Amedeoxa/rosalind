import itertools


def enumeratingGeneOrders(num):
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

if __name__ == '__main__':
    with open('rosalind_perm.txt') as file:
        num = int(file.readline())
    enumeratingGeneOrders(num)