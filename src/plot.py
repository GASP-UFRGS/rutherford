from card_reader import read_card, _raise_missing_card_error 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json
import sys

# Normalizes simulation according to the first datapoint 
def normalize_sim(dataX, dataY, simulationX, simulationY):
	
	minAngleDif = 180 # Arbitrary high initial value

	# Search for simulation point with X closest to smallest value of dataX
	angle_to_find = dataX[0]
	for i in range(len(simulationX)):
		angleDif = abs(simulationX[i] - angle_to_find )
		if angleDif < minAngleDif:
			minAngleDif = angleDif
			minIndex = i
 
	closestAngle = simulationY[minIndex]

	normalizationFactor = dataY[0]/closestAngle
	simulationY = [point*normalizationFactor for point in simulationY]

	return simulationY


# Converts data from nanoBarns to cm^2/sterad to match plots from 1957 Hoftstadt paper
# Converts x axis to cosing if necessary
def convert_data(dataX, dataY, dataYerr, cross_section_variable):

	if cross_section_variable == 'cos':
		dataX = np.cos(np.radians(dataX))

	normalizationFactor = 1e-33

	dataY = [point*normalizationFactor for point in dataY]

	if dataYerr != None:
		dataYerr = [entry * normalizationFactor for entry in dataYerr]

		return dataX, dataY, dataYerr

	return dataX, dataY


def plot(output,card_name):

	#Read parameter input card 
	parameters = read_card(card_name)
	kinEn = parameters.get('kinEn')
	angUnit = parameters.get('angUnit')
	mott = parameters.get('mott')
	recoil = parameters.get('recoil')
	impactParameter = parameters.get('impactParameter') 
	cross_section_variable = parameters.get('difCrossSec')
	hof25 = parameters.get('hoftstadter25') 
	hof125 = parameters.get('hoftstadter125') 
	hof300 = parameters.get('hoftstadter300')
	hof400 = parameters.get('hoftstadter400')
	hof550 = parameters.get('hoftstadter550')
	geiger = parameters.get('GeigerMarsden')
	hof = any([hof25, hof125, hof300, hof400, hof550])
	
	# Read file
	data = pd.read_csv(output)

	# Save each column of output in a dictionary with the header as keys
	column_lists = {}
	for column in data.columns:
	    column_lists[column] = data[column].tolist()

	# Place values on variables
	theta_in = column_lists['theta']

	difCrossSec = {}

	if cross_section_variable in ['cos', 'theta', 'omega']:
		difCrossSec['Rutherford'] = column_lists['difCrossSec_Ruth']
		if mott == "true":
			difCrossSec['Mott'] = column_lists['difCrossSec_Mott']
		if recoil == "true":
			difCrossSec['Target Recoil'] = column_lists['difCrossSec_Recoil']

	if impactParameter == 'true':
		b_out = column_lists['b_out']

		# b vs theta
	
		plt.figure(figsize=(8,6), facecolor='w')
		plt.plot(theta_in, b_out)
		plt.ylabel(r'$b$ [fm]',fontsize=14)
		plt.xlabel(r'$\theta$ [{unit}]'.format(unit=angUnit),fontsize=14)
		plt.title('Impact parameter as function of the scattering angle',fontsize=16)
		pltName = 'plot_b_vs_theta.png'
		plt.savefig(pltName, dpi=300, bbox_inches='tight')
		print(f'Created {pltName}')
	

		# theta vs b

		plt.figure(figsize=(8,6), facecolor='w')
		plt.plot(b_out, theta_in)
		plt.ylabel(r'$\theta$ [{unit}]'.format(unit=angUnit),fontsize=14)
		plt.xlabel(r'$b$ [fm]',fontsize=14)
		plt.title('Scattering angle as function of the impact parameter',fontsize=16)
		pltName = 'plot_theta_vs_b.png'
		plt.savefig(pltName, dpi=300, bbox_inches='tight')
		print(f'Created {pltName}')



	# differential cross section vs dtheta/dcos(theta)/domega

	if cross_section_variable == 'cos':
		theta_in = np.cos(np.radians(theta_in))

	pltName='plot_dsig_'
	plt.figure(figsize=(8,6), facecolor='w')
	plt.yscale("log")


	# Change axis labels
	if cross_section_variable == 'cos':
		pltName += 'dcostheta_vs_costheta'
		plt.xlabel(r'$cos(\theta)$'.format(unit=angUnit),fontsize=14)
		plt.ylabel(r'$d\sigma/dcos(\theta)$',fontsize=14)
		plt.title(r'Distribution of $d\sigma/dcos(\theta)$ as function of the scattering angle',fontsize=16)

	elif cross_section_variable == 'theta': 
		pltName += 'dtheta_vs_theta'
		plt.xlabel(r'$\theta$ [{unit}]'.format(unit=angUnit),fontsize=14)
		plt.ylabel(r'$d\sigma/d\theta$',fontsize=14)
		plt.title(r'Distribution of $d\sigma/d\theta$ as function of the scattering angle',fontsize=16)

	elif cross_section_variable == 'omega':
		pltName += 'domega_vs_theta'
		plt.xlabel(r'$\theta$ [{unit}]'.format(unit=angUnit),fontsize=14)
		plt.ylabel(r'$d\sigma/d\Omega$',fontsize=14)
		plt.title(r'Distribution of $d\sigma/d\Omega$ as function of the scattering angle',fontsize=16)


	# Read Hofstadter data 
	if any([hof300, hof400, hof550]):
		json_file = open(sys.path[0] + '/../data/Hofstadter.json', 'r')
		real = json.load(json_file)
	else:
		#plot every Differential cross section
		for key in difCrossSec:
			plt.plot(theta_in, difCrossSec[key], label=key)
			

	# Extract the x values from the JSON data
	hofAngles = [float(entry["x"][0]["value"]) for entry in real["values"]]

	# Get energy value for labels
	col_names = [E["value"] for E in real["qualifiers"]['E']]


	# Plotting of Hoftstadter data
	if hof25 == 'true':
		if kinEn != 25e6:
			print(f"Chosen kinetic energy {int(kinEn*1e-6)}MeV is not 25 MeV, therefore data was not plotted.")
		else:
			hof_difCrossSec_25 = []
			hofAngles25 = []
			with open(sys.path[0] + '/../data/hoftstadter25.csv', 'r') as hof25:
				for line in hof25:
					angle, value = line.strip().split(';')
					hofAngles25.append(float(angle))
					hof_difCrossSec_25.append(float(value))

			# Converts data units
			hofAngles25, hof_difCrossSec_25= convert_data(hofAngles25, hof_difCrossSec_25, None, cross_section_variable)

			# Normalizing simluation
			for key in difCrossSec:
				difCrossSec[key] = normalize_sim(hofAngles25, hof_difCrossSec_25, theta_in, difCrossSec[key])

			plt.scatter(hofAngles25, hof_difCrossSec_25, label="Hoftstadter 25 Mev") 

	if hof125 == 'true':
		if kinEn != 125e6:
			print(f"Chosen kinetic energy {int(kinEn*1e-6)}MeV is not 125 MeV, therefore data was not plotted.")
		else:
			hof_difCrossSec_125 = []
			hofAngles125 = []
			with open(sys.path[0] + '/../data/hoftstadter125.csv', 'r') as hof125:
				for line in hof125:
					angle, value = line.strip().split(';')
					hofAngles125.append(float(angle))
					hof_difCrossSec_125.append(float(value))

			hofAngles125, hof_difCrossSec_125= convert_data(hofAngles125, hof_difCrossSec_125, None, cross_section_variable)
			
			# Normalizing simluation
			for key in difCrossSec:
				difCrossSec[key] = normalize_sim(hofAngles125, hof_difCrossSec_125, theta_in, difCrossSec[key])
				
			plt.scatter(hofAngles125, hof_difCrossSec_125, marker='s', color='black',  label="Hoftstadter 125 Mev") 


	if hof300 == 'true':
		if kinEn != 300e6:
			print(f"Chosen kinetic energy {int(kinEn*1e-6)}MeV is not 300 MeV, therefore data was not plotted.")
		else:
			hof_difCrossSec_300 = [float(entry["y"][0]["value"]) for entry in real["values"]]
			yerr = [float(entry["y"][2]["errors"][0]["symerror"]) for entry in real["values"]]

			# Converts data units
			hofAngles, hof_difCrossSec_300, yerr = convert_data(hofAngles, hof_difCrossSec_300, yerr, cross_section_variable)

			# Normalizing simluation
			for key in difCrossSec:
				difCrossSec[key] = normalize_sim(hofAngles, hof_difCrossSec_300, theta_in, difCrossSec[key])

			plt.errorbar(hofAngles, hof_difCrossSec_300, yerr, capsize = 3, ls='none', label="Hoftstadter "+col_names[0]) 
		

	if hof400 == 'true':
		if kinEn != 400e6:
			print(f"Chosen kinetic energy {int(kinEn*1e-6)}MeV is not 400 MeV, therefore data was not plotted.")
		else:
			hof_difCrossSec_400 = [float(entry["y"][1]["value"]) for entry in real["values"]]
			yerr = [float(entry["y"][2]["errors"][0]["symerror"]) for entry in real["values"]]

			# Convert data units
			hofAngles, hof_difCrossSec_400, yerr = convert_data(hofAngles, hof_difCrossSec_400, yerr, cross_section_variable)

			# Normalize simluation
			for key in difCrossSec:
				difCrossSec[key] = normalize_sim(hofAngles, hof_difCrossSec_400, theta_in, difCrossSec[key])

			plt.errorbar(hofAngles, hof_difCrossSec_400, yerr, capsize = 3, ls='none', label="Hoftstadter "+col_names[1])
	

	if hof550 == 'true':
		if kinEn != 550e6:
			print(f"Chosen kinetic energy {int(kinEn*1e-6)}MeV is not 550 MeV!, therefore data was not plotted.")
		else:
			hof_difCrossSec_550 = [float(entry["y"][3]["value"]) for entry in real["values"]]
			yerr = [float(entry["y"][3]["errors"][0]["symerror"]) for entry in real["values"]]

			# Convert data units
			hofAngles, hof_difCrossSec_500, yerr = convert_data(hofAngles, hof_difCrossSec_500, yerr, cross_section_variable)

			# Normalize simluation
			for key in difCrossSec:
				difCrossSec[key] = normalize_sim(hofAngles, hof_difCrossSec_500, theta_in, difCrossSec[key])

			plt.errorbar(hofAngles, hof_difCrossSec_550, yerr, capsize = 3, ls='none', label="Hoftstadter "+col_names[3])


	if any([hof25, hof125, hof300, hof400, hof550]):
		pltName += '_hoftstadter'


	if geiger == 'true':
		if kinEn != 125e6:
			print(f"Chosen kinetic energy {int(kinEn*1e-6)}MeV is not 125 MeV!, therefore data was not plotted.")
		else:
			geiger_difCrossSec_125 = []
			geigerAngles125 = []
			with open(sys.path[0] + '/../data/GeigerMarsden.csv', 'r') as hof125:
				for line in hof125:
					angle, value = line.strip().split(';')
					geigerAngles125.append(float(angle))
					geiger_difCrossSec_125.append(float(value))
	
			# Convert data units
			geigerAngles125, geiger_difCrossSec_125 = convert_data(geigerAngles125, geiger_difCrossSec_125, None, cross_section_variable)
			
			# Normalizing simluation
			for key in difCrossSec:
				difCrossSec[key] = normalize_sim(geigerAngles125, geiger_difCrossSec_125, theta_in, difCrossSec[key])
				
			plt.scatter(geigerAngles125, geiger_difCrossSec_125, color='black',  label="GeigerMarsden") 

			pltName += '_GeigerMarsden'


	# Plot every Differential cross section
	for key in difCrossSec:
		plt.plot(theta_in, difCrossSec[key], label=key)
		

	pltName += f'_{int(kinEn*1e-6)}MeV'
	plt.legend()
	plt.savefig(pltName, dpi=300, bbox_inches='tight')
	print(f'Created {pltName}.png')


if __name__ == '__main__':
    plot('output.dat','input.dat')