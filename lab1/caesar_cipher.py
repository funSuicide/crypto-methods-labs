alphabet = 'abcdefghijklmnopqrstuvwxyz'


def caesar_operation(input_file, output_file, key, mode):
    for line in input_file:
        for char in line:
            if char == ' ' or char == '\n':
                output_file.write(char)
            else:
                if mode == 'e':
                    new_index = (alphabet.index(char) + key) % len(alphabet) + len(alphabet) % len(alphabet)
                else:
                    new_index = (alphabet.index(char) - key) % len(alphabet) + len(alphabet) % len(alphabet)
                new_char = alphabet[new_index]
                output_file.write(new_char)
