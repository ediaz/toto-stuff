


def greatestNumber(n, m, k):
    ll = str(n)
    a = [int(i) for i in list(ll)]
    l = a.index(max(a[0:]))
    return ll[l:k-m]


print greatestNumber(97730, 3, 5)
print greatestNumber(4288559, 3, 7)

