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

        '''
        if f_x_sp.is_irreducible:
            return array_a
            '''
        x_sp = sp.Poly.from_list(start_u, gens=x, modulus=p)
        u_x_sp = sp.Poly.from_list(u_x, gens=x, modulus=p)
        for i in range(1, math.floor(m / 2) + 1):
            tmp = sp.Poly.pow(u_x_sp, p)
            #tmp = sp.Poly.from_poly(tmp, modulus=p)
            _, u_x_sp = sp.Poly.div(tmp, f_x_sp)
            #u_x_sp = sp.Poly.from_poly(u_x_sp, modulus=p)
            tmp2 = u_x_sp - x_sp
            d_x = sp.gcd(f_x_sp, tmp2)
            if d_x == 1:
                return array_a

        array_a.clear()


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

    GF256 = galois.GF(7 ** 3)
    print(GF256.properties)

    GF = galois.GF(7 ** 3, repr="poly")
    # test = galois.Poly.Str("x^5 + 143", field=GF)

    f = galois.Poly.Str("2*x^3 + 1", field=GF)
    # f = GF("2*x^3 + 1")
    g = galois.Poly.Str("3*x^3 + 5*x^2 + 2*x + 5", field=GF)
    # g = GF("3*x^3 + 5*x^2 + 2*x + 5")

    print(f*g)

    # s = generate_polynomial(7, 4)
    # print(s)


if __name__ == '__main__':
    main()