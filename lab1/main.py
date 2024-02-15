from caesar_cipher import caesar_operation
from frequency_analysis import frequency_analysis

input_path = 'open_text.txt'
output_path = 'cipher_text.txt'
decrypt_text_path = 'decrypt_text.txt'

open_file = open(input_path, 'r')
cipher_file = open(output_path, 'w')
decrypt_text = open(decrypt_text_path, 'w')

# caesar_operation(open_file, cipher_file, 3, 'e')
frequency_analysis(open_file, decrypt_text)

