from caesar_cipher import caesar_operation
from frequency_analysis import freq_analysis_caesar

plain_text_path = 'alisa.txt'
cipher_text_path = 'cipher_text.txt'
decrypt_text_path = 'decrypt_text.txt'


key = 14


def main():
    plain_text_file = open(plain_text_path, 'r')
    cipher_file = open(cipher_text_path, 'w+')
    decrypt_file = open(decrypt_text_path, 'w')

    caesar_operation(plain_text_file, cipher_file, key, 'e')

    possible_key = freq_analysis_caesar(cipher_file)
    print(possible_key)

    caesar_operation(cipher_file, decrypt_file, possible_key, "d")

    plain_text_file.close()
    cipher_file.close()
    decrypt_file.close()


if __name__ == '__main__':
    main()
