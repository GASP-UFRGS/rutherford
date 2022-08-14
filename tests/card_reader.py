class CardParametersError(Exception):
    '''Exception raised for errors in the parameters of the input card.'''
    pass

# Dictionary that contains the necessary parameters
# of each procedure
procs = {
    'thvsb': ('KinEn', 'Ztarget', 'Zproj')
}

# Empty dictionary that will receive the parameters
# of the procedure to be done
params = {}

def read_card(file_name):
    with open(file_name, 'r') as card:
        for line in card:
            if '[proc]' in line:
                procedure = line.split()[0]
            elif 'parameters:' in line:
                break
        
        for line in card:
            for parameter in procs[procedure]:
                if parameter in line:
                    params[parameter] = float(line.split()[0])

    if tuple(params) != procs[procedure]:
        raise(CardParametersError)

    return tuple(params.values())

if __name__ == '__main__':
    print(read_card('template_card.txt'))

