from caesar_cipher import caesar_operation

input_path = 'open_text.txt'
output_path = 'cipher_text.txt'

input_file = open(input_path, 'r')
output_file = open(output_path, 'w')

caesar_operation(input_file, output_file, 3, 'e')
