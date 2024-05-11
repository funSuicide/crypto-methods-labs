import random
import math
from numpy.polynomial import polynomial as P
import sympy as sp
import numpy as np
import galois

# p = pow(2, 8)
# h = 25


p = 9
h = 3


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
    return math.factorial(par_p) / (math.factorial(par_h) * math.factorial(par_p - par_h))


def get_bin_vector(m: int, par_p: int, par_h: int) -> list:
    current_l = par_h
    M = []
    current_m = m

    for i in range(1, p + 1):
        if current_m >= get_binomial_coefficient(par_p=par_p - i, par_h=current_l):
            M.append(1)
            current_m = current_m - get_binomial_coefficient(par_p=par_p - i, par_h=current_l)
            current_l -= 1
        else:
            M.append(0)

    return M


def get_pi(par_p):
    result = [i for i in range(0, par_p)]
    random.shuffle(result)
    return result


def get_keys(par_p, par_h):
    GF = galois.GF(7 ** 3, repr="poly")
    f = galois.Poly.Str("2*x^3 + 1", field=GF)
    g = galois.Poly.Str("3*x^3 + 5*x^2 + 2*x + 5", field=GF)  # опционально подумать все же над генерацией

    g = GF.primitive_elements[-1]  # надо подумать как выбрать именно наш, ну как варик по списку поискать просто

    test = GF.Range(start=par_p, stop=par_p + par_p)

    Ai = test.log(g)  # вычислить Ai (диск. лог) - проверить

    pi = get_pi(par_p=par_p)
    d = random.randint(0, pow(p, h) - 2)

    Ci = []

    for i in range(0, par_p):
        Ci.append((Ai[pi[i]] + d) % (pow(par_p, par_h) - 1))

    return OpenKey(par_p=par_p, par_h=par_h, par_c=Ci), CloseKey(par_fx=f, par_gx=g, par_pi=pi, par_d=d)


def encrypt(m: int, A: OpenKey):
    en_p = A.p
    en_h = A.h
    en_c = A.c

    len_bin_str = math.floor(get_binomial_coefficient(par_p=en_p, par_h=en_h))

    M = get_bin_vector(m=m, par_p=en_p, par_h=en_h)

    if len(M) != len_bin_str:
        print('proebali')
        return -1, 0

    result_c = []
    for i in range(0, p):
        result_c += (M[i] * en_c[i]) % (pow(en_p, en_h) - 1)

    return 0, result_c


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

    GF = galois.GF(7 ** 3, repr="poly")
    # test = galois.Poly.Str("x^5 + 143", field=GF)

    f = galois.Poly.Str("2*x^3 + 1", field=GF)
    # f = GF("2*x^3 + 1")
    g = galois.Poly.Str("3*x^3 + 5*x^2 + 2*x + 5", field=GF)
    # g = GF("3*x^3 + 5*x^2 + 2*x + 5")

    '''
    Артем вот тут посмотри
    '''
    # print(f * g)
    # test = GF.Range(start=7, stop=7 ** 3)
    # print(test)
    # g = GF.primitive_elements[-1]
    # print(g)
    # print(len(GF.log(test, g)))

    # s = generate_polynomial(7, 4)
    # print(s)

    open_key, close_key = get_keys(par_p=7, par_h=3)
    en_data = encrypt(22, open_key)
    print(en_data)

if __name__ == '__main__':
    main()
