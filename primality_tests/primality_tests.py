import random, math


## requirements

def fast_pow(a, n):
    n_bin = bin(n)
    v = 1
    k = 2  # top
    while True:
        v *= v
        if n_bin[k] == '1':
            v *= a
        k = k + 1
        if k >= len(n_bin):
            break
    return v

def fast_mod(base, expo, mod):
    result = 1
    firstModulus = base % mod
    iteration = expo.bit_length()
    while iteration >= 0:
        result = (result * result) % mod
        if (expo >> iteration) & 1:
            result = (result * firstModulus) % mod
        iteration -= 1
    return result

def verify_polynomial(a):
    c = a.copy()
    i = -1
    while c[i] == 0:
        c.pop()
    return c

def polynomial_add(a, b):
    deg_plus1 = max(len(a), len(b))
    for _ in range(deg_plus1-len(a)):
        a.append(0)
    for _ in range(deg_plus1-len(b)):
        b.append(0)
    c = [a[i]+b[i] for i in range(deg_plus1)]
    return c

def polynomial_sub(a, b):
    deg_plus1 = max(len(a), len(b))
    for _ in range(deg_plus1-len(a)):
        a.append(0)
    for _ in range(deg_plus1-len(b)):
        b.append(0)
    c = [a[i]-b[i] for i in range(deg_plus1)]
    return c

def polynomial_scalar_mul(c,a):
    return [c*a_i for a_i in a]

def polynomial_mul(a, b):
    deg_plus1 = len(a) + len(b) - 1
    c = [0 for _ in range(deg_plus1)]
    for i in range(len(a)):
        for j in range(len(b)):
            c[i+j] += a[i]*b[j]
    return c

def polynomial_scalar_mod(a, m):
    return verify_polynomial([a_i%m for a_i in a])

def polynomial_div(a, b):
    k = len(a)
    m = len(b)
    if k < m:
        return {'q': 1, 'r': a}
    q = []
    r = a[-m:]
    for i in range(k-m):
        q_i = math.floor(r[-1] / b[-1])
        q.insert(0, q_i)
        r = polynomial_sub(r, polynomial_scalar_mul(q_i, b))
        r.insert(0, a[-m-i-1])
        r.pop()
    q_0 = math.floor(r[-1] / b[-1])
    q.insert(0, q_0)
    r = polynomial_sub(r, polynomial_scalar_mul(q_0, b))
    r.pop()
    return {'q': q, 'r': verify_polynomial(r)}
    
def is_exponential(n):
    iteration = math.floor(math.log2(n)) + 1
    for i in range(2, iteration):
        a_float = n**(1/i)
        a_floor = math.floor(a_float)
        a_ceil = math.ceil(a_float)
        if abs(a_floor - a_float) < abs(a_ceil-a_float):
            a = a_floor
        else:
            a = a_ceil
        if fast_pow(a, i) == n:
            return True
    return False

def find_largest_prime_factor(n):
    m = n
    for divisor in range(2, n):
        while m % divisor == 0:
            m = math.floor(m / divisor)
        if m == 1:
            return divisor
    return 1

def identity(a, n, r):
    lhs = [1]
    poly = [a, 1]
    modulo = [-1] + [0 for _ in range(r-1)] + [1]
    bit_len = n.bit_length()
    if n & 1:
        lhs = polynomial_mul(lhs, poly)
    i = 1
    while i < bit_len:
        print("processing multipulation of poly...")
        poly = polynomial_mul(poly, poly)
        print("processing division of poly...")
        poly = polynomial_div(poly, modulo)['r']
        poly = polynomial_scalar_mod(poly, n)
        if (n >> i) & 1:
            print("processing multipulation of lhs...")
            lhs = polynomial_mul(lhs, poly)
            print("processing division of lhs...")
            lhs = polynomial_div(lhs, modulo)['r']
            lhs = polynomial_scalar_mod(lhs, n)
        i += 1

    rhs_deg = n % r
    rhs = [a%n] + [0 for _ in range(rhs_deg-1)] + [1]

    return lhs == rhs

## prime tests

def isPrime(n):
    iteration = math.floor(n**(1/2)+1)
    for i in range(2, iteration):
        if n % i == 0:
            return False
    return True

def fermat_test(n, k=1):
    if n == 2 or n == 3:
        return 'prime'
    for i in range(k):
        a = random.randrange(2, n-1)
        if fast_mod(a, n-1, n) != 1:
            return 'composite'
    return 'probably prime'

def miller_rabin_test(n, k=1):
    if n % 2 == 0: return 'composite'
    s = 1
    while True:
        d = (n-1)/(2**s)
        if d % 2 != 0:
            d = int(d)
            break
        s += 1
    for i in range(k):
        a = random.randrange(2, n-1)
        b = fast_mod(a, d, n)
        if b != 1:
            bool1 = True
        else:
            bool1 = False
        bool2 = True
        for r in range(s):
            b = fast_mod(a, (2**r)*d, n)
            if b == n-1:
                bool2 = False
                break

        if bool1 and bool2:
            return 'composite', a
    return 'probably prime'

def aks_test(n):
    if is_exponential(n):
        return 'composite'
    r = 2
    while r < n:
        if math.gcd(n,r) != 1:
            return 'composite'
        if isPrime(r):
            q = find_largest_prime_factor(r-1)
            if q >= 4*math.sqrt(r)*math.log(n) and \
               fast_mod(n, math.floor((r-1)/q), r) != 1:
                break
        r += 1
    print('r:',r)
    if r == n:
        return 'prime'
    for a in range(1, math.floor(2*math.sqrt(r)*math.log2(n))+1):
        if not identity(a, n, r):
            return 'composite'
    return 'prime'
