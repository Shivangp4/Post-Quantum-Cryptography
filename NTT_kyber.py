
# KYBER Parameters
n = 256
q = 3329
w = 17 #omega


def inv_mod(a, p, mod):
    # Compute modular inverse of a^p
    a_inv = pow(a, mod - 2, mod)
    return pow(a_inv, p, mod)



# Precompute powers of omega, omega inverse
w_inv = pow(w, q-2, q)


powers_of_w = [0]*n
inv_powers_of_w = [0]*n

for i in range(n):
    powers_of_w[i] = pow(w, i, q)
    inv_powers_of_w[i] = pow(w_inv, i, q)





# Cooley-Tukey NTT algorithm in O(nlogn)
def ntt(inp, n = n, factor = 1):

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
    ntt_even = ntt(even, n2, factor * 2)
    ntt_odd = ntt(odd, n2, factor * 2)

    # Combine results
    out = [0] * n

    # ntt(f(i)) = ntt(f_even(i)) + w^i * ntt(f_odd(i))
    # ntt(f(i + n/2)) = ntt(f_even(i)) - w^i * ntt(f_odd(i))
    for i in range(n2):
        out[i] = ntt_even[i] + powers_of_w[(i * factor)] * ntt_odd[i]
        out[i + n2] = ntt_even[i] - powers_of_w[(i * factor)] * ntt_odd[i]
        out[i] = out[i] % q
        out[i + n2] = out[i + n2] % q

    return out

# Cooley-Tukey INTT algorithm in O(nlogn)
def intt(inp, n = n, factor = 1):
    
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
    intt_even = intt(even, n2, factor * 2)
    intt_odd = intt(odd, n2, factor * 2)

    # half = 0.5
    half = pow(2, q-2, q)    

    # Combine results
    out = [0] * n

    # intt(f(i)) =  2^-1 * y^i * (intt(f_even(i)) + w^-i * intt(f_odd(i)))
    # intt(f(i + n/2)) = 2^-1 * (1^(1/4)) * y^i * (intt(f_even(i)) - w^-i * intt(f_odd(i)))
    for i in range(n2):
        out[i] = (half * intt_even[i]) % q + (half * inv_powers_of_w[i * factor] *  intt_odd[i]) % q
        out[i + n2] = (half * intt_even[i]) % q - (half * inv_powers_of_w[i * factor] *  intt_odd[i]) % q
        out[i] = out[i] % q
        out[i + n2] =  out[i + n2] % q

    return out


def multiplication_kyber(a, b):
    

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

        
    n2 = (int) (n // 2)

    a_ntt = ntt_kyber(a)
    b_ntt = ntt_kyber(b)
    

    c_ntt = [0] * n
    

    # Combination step in polynomial multiplication in Kyber
    # c_2i = a_2i * b_2i + w^2i * a_2i+1 * b_2i+1
    # c_2i+1 = a_2i * b_2i+1 + a_2i+1 * b_2i
    for i in range(n2):
        c_ntt[2 * i] = (a_ntt[2 * i] * b_ntt[2 * i] + powers_of_w[2 * i] * a_ntt[2*i+1]*b_ntt[2*i+1]) % q
        c_ntt[2 * i + 1] = (a_ntt[2 * i] * b_ntt[2 * i + 1] + a_ntt[2*i+1]*b_ntt[2*i]) % q    


    # Inverse NTT to get the result
    c = intt_kyber(c_ntt)
    
    print("Result =", c)

    return c





def ntt_kyber(inp):

    n2 = (int) (n // 2)
    even = [0] * n2
    odd = [0] * n2

    for i in range(n2):
        even[i] = inp[2 * i]  # First half of the polynomial
        odd[i] = inp[2 * i + 1]  # Second half of the polynomial


    # Recursively compute NTT of even and odd parts
    # Size of the polynomial is halved, and hence the factor is doubled (as new root of unity becomes w^2)
    even_ntt = ntt(even, n=n2, factor=2)
    odd_ntt = ntt(odd,  n=n2, factor=2)

    out = [0] * (n)
    
    for i in range(n2):
        out[2*i] = (even_ntt[i]) % q
        out[2*i+1] = (odd_ntt[i]) % q
        
    return out


def intt_kyber(inp):


    n2 = (int) (n // 2)
    even = [0] * n2
    odd = [0] * n2

    for i in range(n2):
        even[i] = inp[2 * i]  # First half of the polynomial
        odd[i] = inp[2 * i + 1]  # Second half of the polynomial


    # Recursively compute INTT of even and odd parts
    # Size of the polynomial is halved, and hence the factor is doubled (as new root of unity becomes w^2)
    even_intt = intt(even, n=n2, factor=2)
    odd_intt = intt(odd, n=n2, factor=2)

    out = [0] * (n)

    for i in range(n2):
        out[2*i] = (even_intt[i]) % q
        out[2*i+1] = (odd_intt[i]) % q

    return out




# Test cases

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


multiplication_kyber(a,b)