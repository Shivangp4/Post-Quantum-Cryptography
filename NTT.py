
# Parameters
n = 512
q = 12289
y = 10968  #gamma


def inv_mod(a, p, mod):
    # Compute modular inverse of a^p
    a_inv = pow(a, mod - 2, mod)
    return pow(a_inv, p, mod)



# Precompute powers of omega, omega inverse, gamma, gamma inverse
w = pow(y, 2, q)
w_inv = pow(w, q-2, q)
y_inv = pow(y, q-2, q)

powers_of_w = [0]*n
inv_powers_of_w = [0]*n
powers_of_y = [0]*n*2
inv_powers_of_y = [0]*n*2

for i in range(n):
    powers_of_w[i] = pow(w, i, q)
    inv_powers_of_w[i] = pow(w_inv, i, q)

for i in range(2*n):
    powers_of_y[i] = pow(y, i, q)
    inv_powers_of_y[i] = pow(y_inv, i, q)





# ------------ NAIVE O(n^2) APPROACH ------------

# Naive NTT algorithm in O(n^2)
def ntt(inp):

    out = [0] * n
    for i in range(0, n):
        for j in range(0, n):
            out[i] = out[i] + (inp[j] * powers_of_w[i*j % n] * powers_of_y[ j % (2*n)]) % q 
        out[i] = out[i] % q
    return out

# Naive INTT algorithm in O(n^2)
def intt(inp):
         
    n_inv = pow(n, q-2, q)
    out = [0] * n
    for i in range(0, n):
        for j in range(0, n):
            out[i] = out[i] + (inp[j] * inv_powers_of_w[ i*j % n] * n_inv * inv_powers_of_y[ i % (2*n)] ) % q
        out[i] = out[i] % q
    return out





# ------------ COOLEY-TUKEY O(nlogn) APPROACH ------------


# Cooley-Tukey NTT algorithm in O(nlogn)
def ntt_ct(inp, n = n, factor = 1):

    """ At each recursive step n is halved
    this must mean that (new_omega)^(n/2) = 1
    Therefore at each recursive step omega is doubled
    similarly for gamma
    factor is used to keep track of the powers of omega and gamma """

    # Base case: empty array
    if( n == 1):
        return inp

    n2 = (int) (n // 2)
    odd = [0] * n2
    even = [0] * n2

    for i in range(n2):
        even [ i ] = inp[ 2 * i]
        odd [ i ] = inp[ 2 * i + 1]

    # Recursively compute NTT of even and odd parts
    ntt_even = ntt_ct(even, n2, factor * 2)
    ntt_odd = ntt_ct(odd, n2, factor * 2)

    # Combine results
    out = [0] * n

    # ntt(f(i)) = ntt(f_even(i)) + w^i * ntt(f_odd(i))
    # ntt(f(i + n/2)) = ntt(f_even(i)) - w^i * ntt(f_odd(i))
    for i in range(n2):
        out[i] = ntt_even[i] + powers_of_y[(factor)] * powers_of_w[(i * factor)] * ntt_odd[i]
        out[i + n2] = ntt_even[i] - powers_of_y[(factor)] * powers_of_w[(i * factor)] * ntt_odd[i]
        out[i] = out[i] % q
        out[i + n2] = out[i + n2] % q

    return out



# Cooley-Tukey INTT algorithm in O(nlogn)
def intt_ct(inp, n = n, factor = 1):
    
    """ At each recursive step n is halved
    this must mean that (new_omega)^(n/2) = 1
    Therefore at each recursive step omega is doubled
    similarly for gamma
    factor is used to keep track of the powers of omega and gamma """
 
    # Base case: empty array
    if( n == 1):
        return inp

    n2 = (int) (n // 2)
    odd = [0] * n2
    even = [0] * n2

    for i in range(n2):
        even [ i ] = inp[ 2 * i]
        odd [ i ] = inp[ 2 * i + 1]

    # Recursively compute NTT of even and odd parts
    intt_even = intt_ct(even, n2, factor * 2)
    intt_odd = intt_ct(odd, n2, factor * 2)

    # half = 0.5
    half = pow(2, q-2, q)    

    # Combine results
    out = [0] * n

    # intt(f(i)) =  2^-1 * y^i * (intt(f_even(i)) + w^-i * intt(f_odd(i)))
    # intt(f(i + n/2)) = 2^-1 * (1^(1/4)) * y^i * (intt(f_even(i)) - w^-i * intt(f_odd(i)))
    for i in range(n2):
        out[i] = (half * intt_even[i]) % q + (half * inv_powers_of_w[i * factor] *  intt_odd[i]) % q
        out[i + n2] = (half * intt_even[i]) % q - (half * inv_powers_of_w[i * factor] *  intt_odd[i]) % q
        out[i] = powers_of_y[i * factor] * out[i] % q
        out[i + n2] = 10810 * powers_of_y[(i) * factor] * out[i + n2] % q

    return out



# To calculate polynomial multiplication of a and b
def multiply(a, b):
    
    # Reducing the polynomials in the ring Z_q[x]/(x^n + 1)
    for i in range(n, len(a)):
        a[i - n] -= a[i]
        a[i] = 0

    for i in range(n, len(b)):
        b[i - n] -= b[i]
        b[i] = 0
    

    # Padding zeros to make the length of a and b equal to n
    if len(a) < n:
        a.extend([0] * (n - len(a)))

    if len(b) < n:
        b.extend([0] * (n - len(b)))


    a_ntt = ntt(a[:n])
    b_ntt = ntt(b[:n])
    a_ntt_ct = ntt_ct(a[:n])
    b_ntt_ct = ntt_ct(b[:n])

    # Pointwise multiplication in NTT domain
    c_ntt = [(a_ntt[i] * b_ntt[i]) % q for i in range(n)]
    c_ntt_ct = [(a_ntt_ct[i] * b_ntt_ct[i]) % q for i in range(n)]

    # Inverse NTT to get the result
    c = intt(c_ntt)
    c_ct = intt_ct(c_ntt_ct)

    print("Result =", c)
    print("Result using Cooley-Tukey =", c_ct)
    print((c == c_ct) and a_ntt == a_ntt_ct and b_ntt == b_ntt_ct)  # Check if the results are equal
    

# TEST CASES

a = [1] * 3
b = [1] * 3
## Result = [1, 2, 3, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0,]

# a = [1, 0, 0, 2, 0, 1]
# b = [3, 0, 1, 0, 1]
## Result = [3, 0, 1, 6, 1, 5, 0, 3, 0, 1,]

# a = [1, 2]
# b = [3, 4]
## c = [3, 10, 8]

# a = [1, 2, 3]
# b = [4, 5]
## c = [4, 13, 22, 15]

# a = [2, 1]
# b = [4, 0, 3]
## Result = [8, 4, 6, 3,]

# print(ntt(a))
# print(ntt_ct(a, n))
multiply(a, b)


