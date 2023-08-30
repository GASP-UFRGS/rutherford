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
    'thvsb': ('kinEn', 'zTarget', 'zProj', 'angUnit', 'angStart', 'angEnd', 'mott', 'recoil', 'diracProton', 'formFactor', 'impactParameter', 'difCrossSec', 'hoftstadter25', 'hoftstadter125', 'hoftstadter300', 'hoftstadter400', 'hoftstadter550', 'geigerMarsden')
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


    ####### Tests #######
    if tuple(params) != procs[procedure]:
        _raise_parameters_error()

    # Test for float 
    for variable in ['kinEn', 'angStart', 'angEnd']:
        try:
            variable = float(params[variable])
        except ValueError:
            raise ValueError(f"{params[variable]} is not a valid value for {variable}. It must be a float.")

    # Test for integer
    for variable in ['zTarget', 'zProj']:
        try:
            variable = int(params[variable])
        except ValueError:
            raise ValueError(f"{params[variable]} is not a valid value for {variable}. It must be an integer.")
    
    # Test for booleans
    for variable in ['mott', 'recoil', 'diracProton', 'formFactor', 'impactParameter', 'hoftstadter25', 'hoftstadter125', 'hoftstadter300', 'hoftstadter400', 'hoftstadter550', 'geigerMarsden']:
        if params[variable] not in ['true', 'false']:
            raise ValueError(f"{params[variable]} is not a valid value for {variable}. It must be 'true' or 'false'.")
        else:
            if params[variable] == 'true':
                params[variable] = True
            elif params[variable] == 'false':
                params[variable] = False

    # Test for valid angle unit
    if params['angUnit'] not in ['radians', 'degrees']:
        raise ValueError(f"{params['angUnit']} is not a valid value for angUnit. It must be 'radians' or 'degrees'.")

    # Test for valid choices for difCrossSec
    if params['difCrossSec'] not in ['theta', 'cos', 'omega', 'none']:
        raise ValueError(f"{params['difCrossSec']} is not a valid value for difCrossSec. It must be 'theta', 'cos', 'omega', or 'none'.")

    return params

if __name__ == '__main__':
    print(read_card('input.dat'))
