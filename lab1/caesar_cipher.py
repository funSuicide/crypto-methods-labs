alphabet = 'abcdefghijklmnopqrstuvwxyz'


def caesar_operation(input_file, output_file, key, mode="d"):
    input_file.seek(0)
    output_file.seek(0)

    if mode == "d":
        key = -1 * key

    for line in input_file:
        for char in line:
            if char == ' ' or char == '\n':
                output_file.write(char)
            else:
                new_index = (alphabet.index(char) + key + len(alphabet)) % len(alphabet)
                new_char = alphabet[new_index]
                output_file.write(new_char)
