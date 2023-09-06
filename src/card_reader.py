import platform

class CardError(Exception):
    '''Exception raised for errors related to the input card.'''
    pass

def _raise_parameters_error(params, procedure):
    missing_elements = set(procedure) - set(params)
    raise CardError(f"Parameters missing: {missing_elements}")

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
    'thvsb': ('kinEn', 'zTarget', 'zProj', 'angUnit', 'angStart', 'angEnd'),
    'xsec': ('kinEn', 'zTarget', 'zProj', 'angUnit', 'angStart', 'angEnd', 'mott', 'recoil', 'diracProton', 'formFactor', 'rosenbluth', 'CrossSecVariable', 'hoftstadter25', 'hoftstadter125', 'hoftstadter300', 'hoftstadter400', 'hoftstadter550', 'geigerMarsden'),
    'beam': ('kinEn', 'zTarget', 'zProj', 'angUnit', 'bMin', 'bMax', 'nProj', 'detectorDistance'),
    'bvsd': ('kinEn', 'zTarget', 'zProj', 'angUnit', 'bMin', 'bMax'),
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
                if f'[{parameter}]' in line:
                    try:
                        params[parameter] = float(line.split()[0])
                    except:
                        params[parameter] = line.split()[0]


    ####### Tests #######
    if set(params) != set(procs[procedure]):
        _raise_parameters_error(params, procs[procedure])

    # Test for float 
    floats = ['kinEn', 'angStart', 'angEnd', 'bMin', 'bMax', 'detectorDistance']
    for variable in list(set(params) & set(floats)): # Intersection between params and floats
        try:
            variable = float(params[variable])
        except ValueError:
            raise ValueError(f"{params[variable]} is not a valid value for {variable}. It must be a float.")

    # Test for integer
    integers = ['zTarget', 'zProj', 'nProj']
    for variable in list(set(params) & set(integers)): # Intersection between params and floats
        try:
            params[variable] = int(params[variable])
        except ValueError:
            raise ValueError(f"{params[variable]} is not a valid value for {variable}. It must be an integer.")
    
    # Test for booleans
    booleans = ['mott', 'recoil', 'diracProton', 'formFactor', 'rosenbluth', 'impactParameter', 'hoftstadter25', 'hoftstadter125', 'hoftstadter300', 'hoftstadter400', 'hoftstadter550', 'geigerMarsden']
    for variable in list(set(params) & set(booleans)): # Intersection between params and booleans
        if params[variable] not in ['true', 'false']:
            raise ValueError(f"{params[variable]} is not a valid value for {variable}. It must be 'true' or 'false'.")
        else:
            # Converts string to boolean
            if params[variable] == 'true':
                params[variable] = True
            elif params[variable] == 'false':
                params[variable] = False
                
    # Test for valid angle unit
    units = ['radians', 'degrees']
    if 'angUnit' in procs[procedure]:
        if params['angUnit'] not in ['radians', 'degrees']:
            raise ValueError(f"{params['angUnit']} is not a valid value for angUnit. It must be 'radians' or 'degrees'.")

    # Test for valid choices for difCrossSec
    if 'CrossSecVariable' in params:
       if params['CrossSecVariable'] not in ['theta', 'cos', 'omega']:
            raise ValueError(f"{params['CrossSecVariable']} is not a valid value for CrossSecVariable. It must be 'theta', 'cos' or 'omega'.")

    params['proc'] = procedure
    return params

if __name__ == '__main__':
    print(read_card('input.dat'))
