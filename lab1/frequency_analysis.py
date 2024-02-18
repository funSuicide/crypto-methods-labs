alphabet = 'abcdefghijklmnopqrstuvwxyz'

freq = [0.0817, 0.0149, 0.0278, 0.0425, 0.1270, 0.0223,
        0.0202, 0.0609, 0.0697, 0.0015, 0.0077, 0.0403,
        0.0241, 0.0675, 0.0751, 0.0193, 0.0010, 0.0599,
        0.0633, 0.0906, 0.0276, 0.0098, 0.0236, 0.0015,
        0.0197, 0.0007]

freq_stats = dict(zip(list(alphabet), freq))


def get_freq(input_file):
    input_file.seek(0)
    symbols = dict.fromkeys(list(alphabet), 0)
    counter = 0
    for line in input_file:
        for char in line:
            if not (char == ' ' or char == '\n'):
                symbols[char.lower()] += 1
                counter += 1

    for char in symbols.keys():
        symbols[char] /= counter

    return symbols


def freq_analysis(input_file, output_file):
    symbols = get_freq(input_file)

    sorted_symbols_freq = dict(sorted(symbols.items(), key=lambda item: item[1], reverse=True))
    sorted_freq_stats = dict(sorted(freq_stats.items(), key=lambda item: item[1], reverse=True))

    substitution = dict(zip(sorted_symbols_freq.keys(), sorted_freq_stats.keys()))
    substitution['\n'] = '\n'
    substitution[' '] = ' '

    input_file.seek(0)
    decrypt_text = ''.join([substitution[ch] for line in input_file for ch in line])
    output_file.write(decrypt_text)


def freq_analysis_caesar(crypt_file):
    symbols = get_freq(crypt_file)

    possible_keys = dict.fromkeys(range(len(alphabet)), 0)

    for i in list(alphabet):
        for j in list(alphabet):
            if abs(symbols[i] - freq_stats[j]) < 0.0005:
                possible_keys[(alphabet.index(i) - alphabet.index(j)) % len(alphabet)] += 1

    key = sorted(possible_keys.items(), key=lambda item: item[1], reverse=True)[0][0]
    return key
