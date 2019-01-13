def fibd(n, k=1):  # rosalind_fibd
    ages = [1] + [0]*(k-1)
    for i in range(n-1):
        ages = [sum(ages[1:])] + ages[:-1]
    print(sum(ages))

if __name__ == '__main__':
    fibd(97, 16)