params = {
    'T': '',
    'Z_t': '',
    'z_p': ''
}
def read_card(file_name):
    with open(file_name, 'r') as card:
        for line in card:
            for parameter in params:
                if parameter in line:
                    params[parameter] = float(line.split()[0])

    return params

if __name__ == '__main__':
    print(read_card('template_card.txt'))