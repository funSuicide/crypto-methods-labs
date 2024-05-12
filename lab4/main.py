import random
import math
from numpy.polynomial import polynomial as P
import sympy as sp
import numpy as np
import galois


class OpenKey:
    def __init__(self, par_p, par_h, par_c):
        self.p = par_p
        self.h = par_h
        self.c = par_c


class CloseKey:
    def __init__(self, par_fx, par_gx, par_pi, par_d):
        self.fx = par_fx
        self.gx = par_gx
        self.pi = par_pi
        self.d = par_d


'''
def generate_polynomial(p, m):
    array_a = []
    while True:
        array_a.append(random.randint(1, p - 1))
        for i in range(1, m):
            a_i = random.randint(0, p - 1)
            array_a.append(a_i)
        array_a.append(1)
        array_a.reverse()

        f_x = array_a
        u_x = [1, 0]
        start_u = [1, 0]
        x = sp.Symbol('x')
        f_x_sp = sp.Poly.from_list(f_x, gens=x, modulus=p)

        
if f_x_sp.is_irreducible:
    return array_a
    
x_sp = sp.Poly.from_list(start_u, gens=x, modulus=p)
u_x_sp = sp.Poly.from_list(u_x, gens=x, modulus=p)
for i in range(1, math.floor(m / 2) + 1):
    tmp = sp.Poly.pow(u_x_sp, p)
    # tmp = sp.Poly.from_poly(tmp, modulus=p)
    _, u_x_sp = sp.Poly.div(tmp, f_x_sp)
    # u_x_sp = sp.Poly.from_poly(u_x_sp, modulus=p)
    tmp2 = u_x_sp - x_sp
    d_x = sp.gcd(f_x_sp, tmp2)
    if d_x == 1:
        return array_a

array_a.clear()
'''


def get_binomial_coefficient(par_p: int, par_h: int) -> float:
    if par_p == 0:
        return 0
    if par_h == 0:
        return 1
    return math.factorial(par_p) / (math.factorial(par_h) * math.factorial(par_p - par_h))


def get_bin_vector(m: int, par_p: int, par_h: int) -> list:
    current_l = par_h
    M = [0] * par_p
    current_m = m

    for i in range(1, par_p + 1):
        if current_m >= get_binomial_coefficient(par_p=par_p - i, par_h=current_l):
            M[i - 1] = 1
            current_m -= get_binomial_coefficient(par_p=par_p - i, par_h=current_l)
            current_l -= 1
        else:
            M[i - 1] = 0

    return M


def get_pi(par_p):
    result = [i for i in range(0, par_p)]
    random.shuffle(result)
    return result


def get_keys(par_p, par_h, is_test: bool):
    GF = galois.GF(par_p ** par_h, repr="poly", irreducible_poly="x^4 + 3*x^3 + 5*x^2 + 6*x + 2",
                   primitive_element="3*x^3 + 3*x^2 + 6")
    f = galois.Poly.Str("x^4 + 3*x^3 + 5*x^2 + 6*x + 2", field=GF)
    g_r = galois.Poly.Str("3*x^3 + 3*x^2 + 6", field=GF)  # опционально подумать все же над генерацией

    g = GF("3*x^3 + 3*x^2 + 6")

    test = GF.Range(start=par_p, stop=par_p + par_p)

    Ai = test.log(base=g)

    if is_test:
        pi = [6, 4, 0, 2, 1, 5, 3]
        d = 1702
    else:
        pi = get_pi(par_p=par_p)
        d = random.randint(0, pow(par_p, par_h) - 2)

    Ci = []

    for i in range(0, par_p):
        Ci.append((Ai[pi[i]] + d) % (pow(par_p, par_h) - 1))

    return OpenKey(par_p=par_p, par_h=par_h, par_c=Ci), CloseKey(par_fx=f, par_gx=g_r, par_pi=pi, par_d=d)


def encrypt(m: int, A: OpenKey) -> int:
    en_p = A.p
    en_h = A.h
    en_c = A.c

    M = get_bin_vector(m=m, par_p=en_p, par_h=en_h)

    result_c = 0
    for i in range(0, en_p):
        result_c += (M[i] * en_c[i])

    return result_c % (pow(en_p, en_h) - 1)


def decrypt(c: int, A: OpenKey, B: CloseKey) -> int:
    de_h = A.h
    de_d = B.d
    de_p = A.p
    de_fx = B.fx
    de_gx = B.gx
    de_pi = B.pi

    r = (c - de_h * de_d) % (pow(de_p, de_h) - 1)  # correct

    u_x = np.power(de_gx, r)
    u_x %= de_fx  # correct

    s_x = u_x + de_fx  # correct

    tmp = galois.factors(s_x)

    Ti = []  # correct

    for poly in tmp[0]:
        if len(poly.nonzero_coeffs) != 2:
            Ti.append(0)
        else:
            Ti.append(int(poly.coeffs[1]))

    M = [0] * de_p  # correct

    for t in Ti:
        M[de_pi.index(t)] = 1

    m = 0
    l = de_h

    for i in range(1, de_p + 1):
        if M[i - 1] == 1:
            m = m + get_binomial_coefficient(de_p - i, l)
            l -= 1

    return int(m)


def task(data: int, par_p: int, par_h: int, is_test: bool):
    open_key, close_key = get_keys(par_p=par_p, par_h=par_h, is_test=is_test)
    open_data = data
    print('open data: ', open_data)
    en_data = encrypt(open_data, open_key)
    print('encrypt data: ', en_data)
    dec_data = decrypt(c=en_data, A=open_key, B=close_key)
    print('decrypt data: ', dec_data)


def main():
    # 1. Выбрать конечное поле Fq с характеристикой p, где q = p^h,
    # p ≥ h, и для которого проблема дискретного логарифмирования разрешима.

    # SageCell:
    # p = 7
    # h = 4
    #
    # R = GF(p)['x']
    #
    # irreducible_polynomials = []
    # for p in R.polynomials(h):
    #     if p.is_irreducible():
    #          irreducible_polynomials.append(p)
    #
    # index = randint(0, len(irreducible_polynomials))
    # print(irreducible_polynomials[index])
    #
    # shuffle(irreducible_polynomials)
    # primitive = 0
    # for p in irreducible_polynomials:
    #     if p.is_primitive():
    #         primitive = p
    #         break
    #
    # print(primitive)

    # GF256 = galois.GF(7 ** 3)
    # print(GF256.properties)

    task(22, 7, 4, 1)
    print('-' * 32)
    task(17, 7, 4, 0)
    print('-' * 32)
    task(26, 7, 4, 0)
    print('-' * 32)


if __name__ == '__main__':
    main()
