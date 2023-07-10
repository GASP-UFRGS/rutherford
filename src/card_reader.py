import platform

class CardError(Exception):
    '''Exception raised for errors related to the input card.'''
    pass

def _raise_parameters_error():
    raise CardError("Parameters missing or in wrong order")

def _raise_missing_card_error():
    sys = platform.system()
    if sys in ('Linux', 'Darwin'):
        raise CardError("Missing input card. Try running like "
                        "'python3 core.py input.dat'") from None
    else:
        raise CardError("Missing input card. Try running like "
                        "'python core.py input.dat'") from None

# Dictionary that contains the necessary parameters
# of each procedure
procs = {
    'thvsb': ('kinEn', 'zTarget', 'zProj', 'angUnit', 'angStart', 'angEnd', 'mott')
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
                    try:
                        params[parameter] = float(line.split()[0])
                    except:
                        params[parameter] = line.split()[0]

    if tuple(params) != procs[procedure]:
        _raise_parameters_error()

    return tuple(params.values())

if __name__ == '__main__':
    print(read_card('input.dat'))

