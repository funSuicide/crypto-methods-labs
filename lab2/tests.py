import math
import numpy as np
from scipy.stats import norm


def Maurer_test(data, L, Q, n, expected_value, var_L):
    n = len(data)
    K = int(np.floor(n / L) - Q)

    table = dict.fromkeys([i for i in range(np.power(2, L))], 0)

    for i in range(Q):
        value = int(data[i*L: i*L+L], 2)
        table[value] = i+1

    log_sum = 0
    for i in range(Q, Q+K):
        value = int(data[i*L: i*L+L], 2)
        log_sum += np.log2(i + 1 - table[value])
        table[value] = i + 1

    fn = log_sum / K
    c = 0.7 - 0.8 / L + (4 + 32 / L) * np.power(K, -3 / L) / 15
    sigma = c * np.sqrt(var_L / K)

    p_value = math.erfc(np.abs((fn - expected_value) / (np.sqrt(2)*sigma)))
    print(f"Значение P-value: {p_value}")

    if p_value > 0.01:
        print('Последовательность является случайной')
    else:
        print('Последовательность не является случайной')


def cumulative_sums_test(data, n, mode=0):
    cumsums = np.empty(n)

    cumsums[-1] = 0
    if mode == 1:
        data = data[::-1]

    for i, elem in enumerate(data):
        cumsums[i] = cumsums[i-1] + int(elem) * 2 - 1

    z = np.max(np.abs(cumsums))

    sum1 = 0
    sum2 = 0

    left_border = int((-n/z + 1) / 4)
    right_border = int((n/z - 1) / 4)
    for k in range(left_border, right_border+1):
        sum1 += norm.cdf((4 * k + 1) * z / np.sqrt(n)) - norm.cdf((4*k - 1) * z / np.sqrt(n))

    left_border = int((-n / z - 3) / 4)
    right_border = int((n / z - 1) / 4)
    for k in range(left_border, right_border+1):
        sum2 += norm.cdf((4 * k + 3) * z / np.sqrt(n)) - norm.cdf((4 * k + 1) * z / np.sqrt(n))

    p_value = 1 - sum1 + sum2
    print(f"Значение P-value: {p_value}")

    if p_value > 0.01:
        print('Последовательность является случайной')
    else:
        print('Последовательность не является случайной')


def random_excursions_variant_test(data, n):

    cumsums = np.empty(n+2)
    cumsums[0] = 0
    cumsums[-1] = 0

    for i, elem in enumerate(data):
        cumsums[i+1] = cumsums[i] + int(elem) * 2 - 1

    states = np.arange(1, 10)
    neg_states = -1 * states
    table = dict.fromkeys(np.concatenate((states, neg_states)), 0)

    for key in table.keys():
        table[key] = np.count_nonzero(cumsums == key)

    J = np.count_nonzero(cumsums == 0) - 1

    is_random = True
    for key, value in table.items():
        p_value = math.erfc(np.abs(value - J) / np.sqrt(2 * J * (4 * np.abs(key) - 2)))
        print(f"Значение P-value для {key}: {p_value}")
        if p_value < 0.01:
            is_random = False
            break

    if is_random:
        print('Последовательность является случайной')
    else:
        print('Последовательность не является случайной')


def read_bin(filename):
    with open(filename, "r") as f:
        data = f.read()
    data = data.replace(" ", "")
    data = data.replace("\n", "")
    return data


def main():
    # ≥ 387,840    ≥ 904,960     ≥ 2,068,480
    L = [6, 7, 8]
    Q = [10*2**l for l in L]
    expected_value = [5.2177052, 6.1962507, 7.1836656]
    variance = [2.954, 3.125, 3.238]

    data = read_bin("test.txt")
    n = len(data)
    print(f"Длина последовательности: {n}")
    print("\n\nУниверсальный Маурера")
    Maurer_test(data, L[1], Q[1], n, expected_value[1], variance[1])
    print("\n\nКумулятивных сумм")
    cumulative_sums_test(data, n, 0)
    print("\n\nПроизвольные отклонения (вар. 2)")
    random_excursions_variant_test(data, n)


if __name__ == '__main__':
    main()
