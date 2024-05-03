#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <bitset>
#include <map>

struct openKey {
public:

    int p;
    int h;

    std::vector<int> c;

};


struct closeKey {
public:

    int d;

    std::vector<int> f;
    std::vector<int> g;
    std::vector<int> pi;

};


unsigned long int modexp(unsigned long int x,
    unsigned long int y, unsigned long int N)
{
    if (y == 0) return 1;
    unsigned long int z = modexp(x % N, y / 2, N) % N;
    if (y % 2 == 0)
        return (z * z) % N;
    else
        return ((x % N) * ((z * z) % N)) % N;
}


unsigned long int findPowerModSolution(int a, int b, int p) {
    for (auto i = 1; i <= p; i++) {
        //auto tmp = modexp(b, i, p);
        auto tmp = (b * i) % p;
        if (tmp == a)
            return i;
    }
    return -1;

}


void print(std::vector<int>& f)
{
    if (f.size() == 0) {
        std::cout << "Is empty!" << std::endl;
        return;
    }

    for (auto i = 0; i < f.size() - 1; i++) {
        std::cout << f[i] << "x^" << f.size() - i - 1 << " + ";
    }
    std::cout << f[f.size() - 1] << '\n';
}


void normalize(std::vector<int>& f)
{
    size_t non_zero_index = 0;
    for (size_t i = 0; i < f.size(); ++i) {
        if (f[i] != 0) {
            non_zero_index = i;
            break;
        }
    }
    f.erase(f.begin(), f.begin() + non_zero_index);
}

std::vector<int> add(const std::vector<int>& lhs,
    const std::vector<int>& rhs) {
    int p = 7;
    int lhsDegree = lhs.size();
    int rhsDegree = rhs.size();
    int degree = std::max(lhsDegree, rhsDegree);

    std::vector<int> copy(lhs);
    std::vector<int> copySub(rhs);

    std::vector<int> res(degree);

    //if (lhsDegree < rhsDegree)
    copy.insert(copy.begin(), degree - lhsDegree, 0);

    //if (rhsDegree < lhsDegree)
    copySub.insert(copySub.begin(), degree - rhsDegree, 0);


    for (auto i = 0; i < lhsDegree; i++) {
        int tmp = (copy[i] + copySub[i]) % p;
        res[i] = tmp;
    }

    normalize(res);
    return res;
}


std::vector<int> sub(const std::vector<int>& lhs,
    const std::vector<int>& rhs) {
    int p = 7;
    int lhsDegree = lhs.size();
    int rhsDegree = rhs.size();
    int degree = std::max(lhsDegree, rhsDegree);

    std::vector<int> copy(lhs);
    std::vector<int> copySub(rhs);

    std::vector<int> res(degree);

    //if (lhsDegree < rhsDegree)
    copy.insert(copy.begin(), degree - lhsDegree, 0);

    //if (rhsDegree < lhsDegree)
    copySub.insert(copySub.begin(), degree - rhsDegree, 0);

    for (auto i = 0; i < lhsDegree; i++) {
        int tmp = copy[i] - copySub[i];
        if (tmp < 0) tmp += p;
        res[i] = tmp;
    }

    normalize(res);
    return res;
}


std::vector<int> mul(const std::vector<int>& lhs,
    const std::vector<int>& rhs) {
    std::vector<int> res(lhs.size() + rhs.size() - 1, 0);
    for (auto i = 0; i < lhs.size(); i++)
        for (auto j = 0; j < rhs.size(); j++)
            res[i + j] += lhs[i] * rhs[j];

    normalize(res);
    return res;
}


std::pair<std::vector<int>, std::vector<int>> div(
    const std::vector<int>& lhs,
    const std::vector<int>& rhs) {

    int p = 7;

    if (lhs.size() < rhs.size()) {
        std::pair<std::vector<int>, std::vector<int>> res;
        std::vector<int> div(0);
        std::vector<int> rem(lhs);
        res.first = div;
        res.second = rem;
        return res;
    }



    std::vector<int> q(lhs);
    std::vector<int> r(0);

    std::map<int, int> result;

    while (q.size() >= rhs.size()) {
        int diffDegree = q.size() - 1 - (rhs.size() - 1);
        int coefficient = findPowerModSolution(q[0], rhs[0], p);

        std::vector<int> subtracted(rhs);

        for (auto i = 0; i < subtracted.size(); i++) {
            subtracted[i] = (coefficient * subtracted[i]) % p;
        }

        //std::transform(subtracted.begin(), subtracted.end(), subtracted.begin(), [coefficient, p](int xx) { return (coefficient * xx) % p; });
        subtracted.insert(subtracted.end(), diffDegree, 0);

        r = sub(q, subtracted);

        result[diffDegree] = coefficient;

        q = r;
    }


    auto maxKey = *std::max_element(result.begin(), result.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
        });

    std::vector<int> res(maxKey.first + 1, 0);
    for (const auto& element : result) {
        res[element.first] = element.second;
    }

    std::reverse(res.begin(), res.end());
    std::pair<std::vector<int>, std::vector<int>> answer(res, r);

    return answer;
}


unsigned long int solve(const std::vector<int>& f, int x) {
    std::vector<int> e = f;
    std::reverse(e.begin(), e.end());

    unsigned long int result = 0.0;
    int power = 1;

    for (int coefficient : e) {
        result += coefficient * power;
        power *= x;
    }

    return result;
}


unsigned long int factorial(int n) {
    return (n == 0) || (n == 1) ? 1 : n * factorial(n - 1);
}


int discreteLogarithm(std::vector<int> g, std::vector<int> f, int p, int q, int coef = 0) {
    std::vector<int> t(g);

    for (auto x = 2; true; x++) {
        t = mul(t, g);

        std::transform(t.begin(), t.end(), t.begin(), [p](int xx) { return  xx % p; });

        std::pair<std::vector<int>, std::vector<int>> d;
        d = div(t, f);
        t = d.second;

        auto tmp = solve(t, x);

        if (tmp == x + coef) {
            return x;
        }

    }
    return 0;
}



/// <summary>
/// 2.218
/// </summary>
/// <param name="g"></param>
/// <param name="h"></param>
/// <returns></returns>
std::vector<int> euclidean(std::vector<int>& g, std::vector<int>& h) {

    std::vector<int> gg(g);
    std::vector<int> hh(h);

    std::vector<int> r;

    while (hh.size() != 0) {

        auto d = div(gg, hh);
        r = d.second;

        gg = hh;
        hh = r;
    }

    return gg;
}


/// <summary>
/// 4.69
/// </summary>
/// <param name="f"></param>
/// <param name="p"></param>
/// <param name="h"></param>
/// <returns></returns>
bool isIrreducible(std::vector<int>& f, int p, int h) {

    std::vector<int> u;

    for (auto i = 1; i < h / 2; i++) {
        // 2.227 заменяю на мою прекрасную функцию умножения

        // что понимается под u(x) ?



    }
    return true;
}


/// <summary>
/// 4.70
/// </summary>
/// <param name="p"></param>
/// <param name="h"></param>
/// <returns></returns>
std::vector<int> monicIrreduciblePolynomial(const int p, const int h) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, p - 1);

    std::vector<int> irreducible(h, 1);

    do {

        for (auto i = 1; i < h; i++)
            irreducible[i] = dis(gen);

        while (irreducible[h - 1] == 0)
            irreducible[h - 1] = dis(gen);

    } while (!isIrreducible(irreducible, p, h));


    std::vector<int> f{ 1, 3, 5, 6, 2 };
    return f;
    //return irreducible;
}


int gcd(int a, int b) {
    while (b != 0) {
        int t = b;
        b = a % b;
        a = t;
    }
    return a;
}
std::vector<int> factorize(int a) {
    std::vector<int> factors; // Вектор для хранения факторов
    for (int i = 2; i <= sqrt(a); ++i) {
        // Пока i меньше квадратного корня из a
        while (a % i == 0) {
            // Если i является делителем a, добавляем его в факторы и делим a на i
            factors.push_back(i);
            a /= i;
        }
    }
    // Если после деления a больше 1, то a - простое число
    if (a > 1) {
        factors.push_back(a);
    }
    return factors;
}


/// <summary>
/// 4.80
/// </summary>
/// <param name="p"></param>
/// <returns></returns>
std::vector<int> findingGeneratorOfACyclicGroup(int p) {
    std::vector<int> fact;

    fact = factorize(p);


    std::vector<int> g{ 3, 3, 0, 6 };
    return g;

}



std::vector<int> permutation(int size) {

    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<int> p(size);
    std::iota(p.begin(), p.end(), 0);
    std::shuffle(p.begin(), p.end(), gen);

    std::vector<int> pp{ 6, 4, 0, 2, 1, 5, 3 };

    return pp;
}


int getD(int begin, int end) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(begin, end);

    return 1702;
    //return dis(gen);
}


std::pair<openKey, closeKey> keyGeneration(int p, int h) {

    int q = std::pow(p, h);

    std::vector<int> f;
    std::vector<int> g;

    f = monicIrreduciblePolynomial(p, h);
    g = findingGeneratorOfACyclicGroup(p);

    std::vector<int> a;

    for (auto x = 0; x < p; x++)
        a.push_back(discreteLogarithm(g, f, p, q, x));

    std::vector<int> pi;
    pi = permutation(p);

    int d = getD(0, q - 2);

    std::vector<int> c(p);
    for (auto i = 0; i < p; i++)
        c[i] = (a[pi[i]] + d) % (q - 1);

    openKey open;
    open.c = c;
    open.h = h;
    open.p = p;

    closeKey close;
    close.d = d;
    close.f = f;
    close.g = g;
    close.pi = pi;

    std::pair<openKey, closeKey> keys(open, close);
    return keys;
}


int binomialCoefficient(int p, int h) {
    return factorial(p) / (factorial(h) * factorial(p - h));
}


std::vector<int> getM(int m, int p, int h) {

    std::vector<int> M(p);
    int l = h;

    for (auto i = 1; i <= p; i++) {

        int tmp = 0;
        if (p != i && l == 0) {
            tmp = 1;
        }
        if (p > i) {
            tmp = binomialCoefficient(p - i, l);
        }

        if (m >= tmp) {
            M[i - 1] = 1;
            m -= tmp;
            l -= 1;
        }
        else {
            M[i - 1] = 0;
        }

        if (l == 0 && tmp >= 0)
            M[i - 1] = 1;
        if (tmp == 0 && l >= 1)
            M[i - 1] = 0;

    }
    return M;
}


int encryption(openKey key, int m) {

    int binomial_coefficient = log2(binomialCoefficient(key.p, key.h));

    auto M = getM(m, key.p, key.h);

    int q = std::pow(key.p, key.h) - 1;

    int c = 0;
    for (auto i = 0; i < key.p; i++)
        c = (c + M[i] * key.c[i]) % q;

    return c;
}


int decryption(std::pair<openKey, closeKey> keys, int c) {

    auto open = keys.first;
    auto close = keys.second;

    int p = open.p;

    int q = std::pow(open.p, open.h) - 1;
    auto r = (c - open.h * close.d) % q;
    if (r < 0) r += q;

    std::vector<int> u(close.g);

    for (auto i = 1; i < r; i++) {
        u = mul(u, close.g);
        std::transform(u.begin(), u.end(), u.begin(), [p](int xx) { return  xx % p; });
        auto tmp = div(u, close.f);
        u = tmp.second;
    }

    std::vector<int> s;
    s = add(u, close.f);

    print(s);
    std::cout << std::endl;


    int size = 0;
    std::cout << "Input solve" << std::endl;
    std::cin >> size;

    std::vector<int> t(size);
    for (auto i = 0; i < size; i++)
        std::cin >> t[i];

    std::vector<int> M(p, 0);
    for (auto i = 0; i < size; i++) {
        M[t[i]] = 1;
    }

    int m = 0;
    int l = open.h;
    for (auto i = 1; i < p; i++) {
        if (M[i - 1] == 1) {
            m += binomialCoefficient(p - i, l);
            l--;
        }

    }

    return m;
}




int main()
{
    int p = 7;
    int h = 4;

    auto keys = keyGeneration(p, h);

    auto open = keys.first;
    auto close = keys.second;

    auto encryptionData = encryption(open, 22);

    std::cout << "Encryption data: " << encryptionData << std::endl;

    auto decryptionData = decryption(keys, encryptionData);

    std::cout << "\nDecryption data: " << decryptionData << std::endl;

    std::cout << "\n\nHello world :)";


}

/*
    std::bitset<5> bs(m);

    std::vector<int> bits;
    for (auto i = 0; i < bs.size(); i++)
        bits.push_back(bs[i]);

    std::reverse(bits.begin(), bits.end());
*/

/*
    for (auto x = 2; x < q; x++) {
        t = mul(t, g);

        std::transform(t.begin(), t.end(), t.begin(), [p](int xx) { return  xx % p; });

        std::pair<std::vector<int>, std::vector<int>> d;
        d = div(t, f);
        t = d.second;

        auto tmp = solve(t, x);


        if (tmp == x + 3) {
            std::cout << x;
        }



    }

    int p = 7;
    int h = 4;
    int d = 1702;

    std::vector<int> permutation = { 6, 4, 0, 2, 1, 5, 3 };
    std::vector<int> c = { 1925, 2081, 330, 1356, 1237, 1082, 310 };

    std::vector<int> g{ 3, 3, 0, 6 };
    std::vector<int> f{ 1, 3, 5, 6, 2 };


    int binomial_coefficient = log2(factorial(p) / (factorial(h) * factorial(p - h)));

    int message = 22;

    std::bitset<5> bs(message);

    std::vector<int> M{ 1, 0, 1, 1, 0, 0, 1 };

    int cc = (c[0] + c[2] + c[3] + c[6]) % 2400;

    std::cout << cc << std::endl;

    // дешифрование

    int r = (cc - h * d) % 2400;
    if (r < 0) r += 2400;

    std::cout << r << std::endl;
*/



/*
    std::cout << modexp(0, 6, 7) << std::endl;
    std::cout << modexp(3, 6, 7) << std::endl;
    std::cout << modexp(6, 6, 7);

#include "polynomial.h"


//Удаление старших одночленов с нулевыми коэффициентами
void normalize(std::vector<int>& f)
{
    while (f.size() > 1 && f.back() == 0)
        f.pop_back();
}
//Умножение многочленов F(X) * G(X)
std::vector<int> multiple(std::vector<int>& f, std::vector<int>& g)
{
    std::vector<int> h(f.size() + g.size() - 1, 0);
    for (int i = 0; i < f.size(); i++)
        for (int j = 0; j < g.size(); j++)
            h[i + j] += f[i] * g[j];

    normalize(h);
    return h;
}

//Композиция многочленов F(G(X))
std::vector<int> composition(std::vector<int>& f, std::vector<int>& g)
{
    std::vector<int> h = { f.back() };
    for (int i = f.size() - 2; i >= 0; i--)
    {
        h = multiple(h, g);
        h[0] += f[i];
    }

    normalize(h);
    return h;
}

//Хоть сколько-то наглядный вывод многочлена
void print(std::vector<int>& f)
{
    for (int i = f.size() - 1; i > 0; i--)
        std::cout << f[i] << "*X^" << i << " + ";
    std::cout << f[0] << '\n';
}

unsigned long int evaluatePolynomial(const std::vector<int>& coefficients, unsigned long int x) {

    std::vector<int> e = coefficients;
    std::reverse(e.begin(), e.end());

    unsigned long int result = 0.0;
    int power = 1;

    for (int coefficient : e) {
        result += coefficient * power;
        power *= x;
    }

    return result;
}

void permutationLala() {

    const int p = 5; // Размер массива
    std::vector<int> vector(p);
    std::iota(vector.begin(), vector.end(), 0);

    std::random_device rd;
    std::mt19937 g(rd());

    std::shuffle(vector.begin(), vector.end(), g);

    for (auto i = 0; i < p; i++)
        std::cout << vector[i] << std::endl;
}

unsigned long int factorial(int n) {
    return (n == 0) || (n == 1) ? 1 : n * factorial(n - 1);
}

unsigned long int g(unsigned long int x) {
    return 3 * x * x * x + 3 * x * x + 6;
}

unsigned long int f(unsigned long int x) {
    return x * x * x * x + 3 * x * x * x + 5 * x * x + 6 * x + 2;
}

unsigned long int modexp(unsigned long int x, unsigned long int y, unsigned long int N)
{
    if (y == 0) return 1;
    unsigned long int z = modexp(x % N, y / 2, N) % N;
    if (y % 2 == 0)
        return (z * z) % N;
    else
        return ((x % N) * ((z * z) % N)) % N;
}



//std::vector<int> Remainder(std::vector<int>& p1, std::vector<int>& p2)
//{
//    std::vector<int> r;
//    std::vector<int> q;		//Remainder, Quotient
//
//    auto deg = p1.size() - p2.size() - 2;
//    auto rem = std::vector<int>(deg + 1, 0);
//
//    auto p1Deg = p1.size() - 1;
//    auto p2Deg = p2.size() - 1;
//
//    for (int i = 0; p2Deg <= p1Deg; i++)
//    {
//        q[deg - i] = p1[p1Deg] / p2[p2Deg];
//
//        std::vector<int> temp1;
//        auto tmpDeg = deg - i;
//        std::vector<int> temp1(tmpDeg + 1);
//        temp1[deg - i] = q[deg - i];
//
//        Polynomial temp2 = multiply(temp1, p2);
//        r = subtract(p1, temp2);
//        p1 = r;
//    }
//    return r;
//}

//std::vector<int> dividePolynomials(std::vector<int>& dividend, std::vector<int>& divisor) {
//    std::vector<int> quotient;
//
//    int m = dividend.size();
//    int n = divisor.size();
//
//    if (n > m) {
//        quotient.push_back(0);
//        return quotient;
//    }
//
//    std::vector<int> tempDividend = dividend;
//
//    while (m >= n) {
//        int coeff = tempDividend[m - 1] / divisor[n - 1];
//        quotient.insert(quotient.begin(), coeff);
//
//        for (int i = 0; i < n; ++i) {
//            tempDividend[m - n + i] -= coeff * divisor[i];
//        }
//
//        while (m > 0 && tempDividend[m - 1] == 0) {
//            tempDividend.pop_back();
//            m--;
//        }
//    }
//
//    return quotient;
//}



*/

/*
struct Polynomial
{
    int deg;
    float* c;
};

float getCoff(Polynomial p, int term)
{
    if (term > p.deg)
        return 0;
    return p.c[term];
}
void copy(float* a, float* b, int size)
{
    for (int i = 0; i < size; i++)
        a[i] = b[i];
}
void ShrinkPolynomial(Polynomial& p) {

    for (int i = p.deg; i >= 0; i--)
    {
        if (p.c[i] == 0)
            p.deg--;
        else break;
    }
    float* temp = new float[p.deg + 1];

    copy(temp, p.c, p.deg + 1);

    delete[] p.c;
    p.c = temp;
}
Polynomial subtract(Polynomial p1, Polynomial p2)
{
    Polynomial S;
    S.deg = p1.deg >= p2.deg ? p1.deg : p2.deg;
    S.c = new float[S.deg + 1];

    for (int term = 0; term < S.deg + 1; term++)
        S.c[term] = getCoff(p1, term) - getCoff(p2, term);

    ShrinkPolynomial(S);

    return S;
}
Polynomial multiply(Polynomial p1, Polynomial p2)
{
    Polynomial mul;
    mul.deg = p1.deg + p2.deg;
    mul.c = new float[mul.deg + 1] {};

    for (int i = 0; i < p1.deg + 1; i++)
        for (int j = 0; j < p2.deg + 1; j++)
            mul.c[i + j] = mul.c[i + j] + (getCoff(p1, i) * getCoff(p2, j));

    return mul;
}
Polynomial Remainder(Polynomial p1, Polynomial p2)
{
    Polynomial r, q;		//Remainder, Quotient

    r.c = nullptr;
    r.deg = 0;
    q.deg = p1.deg - p2.deg;
    q.c = new float[q.deg + 1] {};

    for (int i = 0; p2.deg <= p1.deg; i++)
    {
        q.c[q.deg - i] = p1.c[p1.deg] / p2.c[p2.deg];

        Polynomial temp1;
        temp1.deg = q.deg - i;
        temp1.c = new float[temp1.deg + 1] {};
        temp1.c[q.deg - i] = q.c[q.deg - i];

        Polynomial temp2 = multiply(temp1, p2);
        r = subtract(p1, temp2);
        p1 = r;
    }
    return r;
}
*/

/*


    int p = 7;
    int h = 4;
    int m = 1702;

    std::vector<int> permutation = {6, 4, 0, 2, 1, 5, 3};

    std::vector<int> c = { 1925, 2081, 330, 1356, 1237, 1082, 310 };

    // C(p, h) = p! / (h!  *  (p - h)!)
    int b = log10(factorial(p) / (factorial(h) * factorial(p - h)));

    std::vector<int> polynomial11{ 6, 0, 3, 3 };
    std::vector<int> polynomial22{ 2, 6, 5, 3, 1 };


    std::vector<int> quotient = dividePolynomials(polynomial11, polynomial22);

    // Вычисляем остаток
    std::vector<int> remainder = polynomial11;
    for (int coeff : quotient) {
        for (int i = 0; i < remainder.size(); ++i) {
            remainder[i] -= coeff * polynomial22[i];
        }
    }

    std::cout << "Remainder: ";
    for (int coeff : remainder) {
        std::cout << coeff << " ";
    }


    for (auto i = 1; i < m; i++) {

        //polynomial11 = multiple(polynomial11, polynomial11);

        auto ff = evaluatePolynomial(polynomial22, i);
        auto gg = evaluatePolynomial(polynomial11, i);

        auto gx = modexp(gg, i, ff);

        if (gx == i) {
            std::cout << gx << "\t\t" << i << std::endl;
            break;
        }

    }
*/


/*

    for (auto x = 0; x < 1; x++) {
        auto ff = f(x);
        polynomial1 = multiplyPolynomials(polynomial1, polynomial1);

        std::cout << std::endl;
        for (auto coeff : polynomial1) {
            std::cout << coeff << " ";
        }
        std::cout << std::endl;

        auto gg = evaluatePolynomial(polynomial1, x);


        auto gx = modexp(gg, x, ff);

        std::cout << x << "\t\t" << gg << "\t\t" << gx << "\t\t" << f << std::endl;

        if (gx == x) {
            std::cout << gx << "\t\t" << x << std::endl;
            break;
        }


    }

*/

/*
/// <summary>
/// example
/// m = 4: x^4 + a1 x^3 + a2 x^2 + a3 x^1 + a4
/// in vector: x[4] = 1, x[3] = a1, x[2] = a2, x[1] = a3, x[0] = a4
/// </summary>
class polynomial {
private:

public:

    /// <summary>
    /// maximum degree of the polynomial
    /// </summary>
    int _degree = 0;

    /// <summary>
    /// container storing pairs of type (int, int)
    /// where the first value is an element of the polynomial,
    /// the second is its value
    /// </summary>
    std::vector<int> _coefficients;

    polynomial() : _coefficients() {}

    polynomial(const polynomial& other) {
        _degree = other._degree;
        _coefficients = other._coefficients;
    }

    polynomial(const int degree) {
        _degree = degree;
        _coefficients = std::vector<int>(_degree + 1, 0);
    }


    /// <summary>
    /// 4.69
    /// </summary>
    /// <param name="prime"></param>
    /// <returns></returns>
    bool testPolynomialForIrreducible(int prime = 7) {

    }

    /// <summary>
    /// 4.70
    /// </summary>
    /// <param name="prime"></param>
    void randomMonicIrreduciblePolynomial(int prime = 7) {
        // randomly select integers between 0 and p-1, with a0 != 0
        std::generate(
            _coefficients.begin() + 1,
            _coefficients.end() - 1,
            std::rand() % prime);

        _coefficients[_degree] = 1;

        do {
            _coefficients[0] = std::rand() % prime;
        } while (_coefficients[0] == 0);

    }


};


int repeatedSquareAndMultiply() {

}





/// <summary>
/// Testing a polynomial for irreducibility
/// </summary>
/// <param name="p">prime number</param>
/// <param name="polynomial"></param>
/// <returns></returns>
bool isIrreducible(int p, polynomial polynomial) {


    return true;
}

/// <summary>
/// 4.70. Generating a random monic irreducible polynomial over Zp
/// </summary>
/// <param name="p">prime number</param>
/// <param name="m">positive integer</param>
/// <returns></returns>
polynomial randomIrreduciblePolynomial(int p, unsigned m) {


}




void keyGeneration() {

    // 1. Выбрать конечное поле Fq с характеристикой p, где q = p^h,
    // p ≥ h, и для которого проблема дискретного логарифмирования разрешима.
    int p = 7, h = 4;

    // 2. Выбрать случайный неприводимый многочлен f(x) степени h над Zp
    // (см.Приложение п.4.70).Элементы Fq должны быть представлены в виде
    // полиномов в Zp[x] степени, меньшей h, с операцией умножения по модулю f(x).

    polynomial f = polynomial(4);
    f._coefficients[0] = 2;
    f._coefficients[1] = 6;
    f._coefficients[2] = 5;
    f._coefficients[3] = 3;
    f._coefficients[4] = 1;

    // 3. Выбрать случайный примитивный многочлен g(x) поля Fq (см. Приложение п.4.80).

    polynomial g = polynomial(3);
    f._coefficients[0] = 6;
    f._coefficients[1] = 0;
    f._coefficients[2] = 3;
    f._coefficients[3] = 3;

    // 4. Для каждого элемента поля i ∈ Zp, найти дискретный логарифм ai = log_g(x) (x+i).

    // 5. Выбрать случайную перестановку π множества {0, 1, 2,... ,p − 1}.

    // 6. Выбрать случайное целое d, 0 ≤ d ≤ p^h − 2.

    // 7. Вычислить ci = (a__π(i) + d) mod (p^h − 1), 0 ≤ i ≤ p − 1.

    // 8. ((c0, c1 ,... , cp−1), p, h) – открытый ключ; (f(x), g(x), π, d) – закрытый ключ.

}

*/