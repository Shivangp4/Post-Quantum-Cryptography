# Post-Quantum Cryptography

As a part of CS674 post-quantum cryptography, I implemented the Number Theoretic Transform (NTT) according to NewHope protocol using Cooley-Tukey divide and conquer approach. NTT (and inverse NTT) was used to calculate polynomial multiplication in Ring $R_q/x^n+1$ in $O(n\log n)$ time. Also, implemented NTT and polynomial multiplication in Crsytal-Kyber protocol.


## Regular approach

NTT is defined as

$$ \hat{f_i} = \sum_{j} f_j \omega^{ij}$$

where { $\hat{f_i}$ } is the NTT of the polynomial { $f_j$ }, and $\omega$ is the n-th root of unity, where n is the degree of polynomial $f$
Inverse NTT (INTT) is defined as 

$$ f_i = \sum_{j} \hat{f_j} \omega^{-ij}$$

## NewHope approach

$$ \hat{f_i} = \sum_{j} f_j \gamma^{j} \omega^{ij}$$
$$ f_i = n^{-1} \gamma^{-i} \sum_{j} \hat{f_j} \omega^{-ij}$$

Polynomial multiplication is defined by $h(x) = f(x)g(x)$ as point-wise multiplication in NTT domain

$$ \hat{h_i} = \hat{f_i} \times \hat{g_i}$$


## Kyber approach

$$ f(x) = f_{even} (x^2) + xf\{odd}(x^2) $$
$$ NTT(f(x)) = NTT(f_{even}(x)) + x NTT(f_{odd}(x)) $$
$$ INTT(f(x)) = INTT(f_{even}(x)) + x INTT(f_{odd}(x)) $$

Polynomial multiplication is defined as 

$$ \hat{h_{2i}} = \hat{f_{2i}} \hat{g_{2i}} + \omega^{2i} \hat{f_{2i+1}} \hat{g_{2i+1}} $$

$$ \hat{h_{2i+1}} = \hat{f_{2i+1}} \hat{g_{2i}} + \hat{f_{2i}} \hat{g_{2i+1}} $$ 

Newhope (and Cooley-Tukey approach) is implemented in """ntt.py""" and kyber approach is implemented in """ntt_kyber.py"""
