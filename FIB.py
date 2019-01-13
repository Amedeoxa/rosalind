def rabbitsAndReacurrence(n, k=1): ## fibonacci
    memo = {0: 0, 1: 1}
    if n not in memo:
        memo[n] = (rabbitsAndReacurrence(n - 1, k) + rabbitsAndReacurrence(n - 2, k)*k)
    return memo[n]