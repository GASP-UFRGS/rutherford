from collections import OrderedDict

params = {
    'T_c': '',
    'Z_t': '',
    'z_p': ''
}

def read_card(file_name):
    with open(file_name, 'r') as card:
        for line in card:
            for parameter in params:
                if parameter in line:
                    params[parameter] = float(line.split()[0])

    return tuple(OrderedDict(params).values())

if __name__ == '__main__':
    print(read_card('template_card.txt'))