alphabet = 'abcdefghijklmnopqrstuvwxyz'
ignore = [' ', ', ', '\n', ';', ':', '\t', '"', '-', '?', '!', '(', ')']



def caesar_operation(input_file, output_file, key, mode="d"):
    input_file.seek(0)
    output_file.seek(0)

    if mode == "d":
        key = -1 * key

    for line in input_file:
        for char in line:
            if char.lower() in alphabet:
                new_index = (alphabet.index(char.lower()) + key + len(alphabet)) % len(alphabet)
                new_char = alphabet[new_index]
                output_file.write(new_char)
            else:
                output_file.write(char)
