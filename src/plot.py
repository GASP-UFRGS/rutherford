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
# Converts x axis to cosine if necessary
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
	procedure = parameters.get('proc')
	kinEn = parameters.get('kinEn')
	angUnit = parameters.get('angUnit')
	bMin = parameters.get('bMin')
	bMax = parameters.get('bMax')
	nProj = parameters.get('nProj')
	detectorDistance = parameters.get('detectorDistance')
	mott = parameters.get('mott')
	recoil = parameters.get('recoil')
	diracProton = parameters.get('diracProton')
	formFactor = parameters.get('formFactor')
	rosenbluth = parameters.get('rosenbluth')
	cross_section_variable = parameters.get('CrossSecVariable')
	hof25 = parameters.get('hoftstadter25') 
	hof125 = parameters.get('hoftstadter125') 
	hof300 = parameters.get('hoftstadter300')
	hof400 = parameters.get('hoftstadter400')
	hof550 = parameters.get('hoftstadter550')
	geiger = parameters.get('geigerMarsden')
	hof = any([hof25, hof125, hof300, hof400, hof550])
	
	# Read file
	data = pd.read_csv(output)

	# Save each column of output in a dictionary with the header as keys
	column_lists = {}
	for column in data.columns:
	    column_lists[column] = data[column].tolist()

	if procedure == 'beam':
		theta_in = column_lists['theta']
		phi = column_lists['phi']
		position = column_lists['position']
		
		# Number of bins
		num_bins = 50

		# Calculate the minimum and maximum values of theta_in
		min_theta = min(theta_in)
		max_theta = max(theta_in)

		# Compute the range and bin width
		theta_range = max_theta - min_theta
		bin_width = theta_range / num_bins

		# Initialize bins and group phi values
		bins = np.linspace(min_theta, max_theta, num_bins+1)
		phi_binned = [0] * num_bins

		# Group phi values into bins
		for i in range(len(theta_in)):
			bin_index = int((theta_in[i] - min_theta) // bin_width)
			bin_index = min(bin_index, num_bins-1)  # Ensure index is within bounds
			phi_binned[bin_index] += phi[i]

		# Calculate bin centers for plotting
		bin_centers = [(bins[i] + bins[i+1]) / 2 for i in range(num_bins)]

		# Create the bar plot
		plt.figure(figsize=(8, 6))
		plt.bar(bin_centers, phi_binned, width=bin_width, align='center')
		plt.xlim(min_theta, max_theta)
		plt.ylabel(r'$\phi$ [{unit}]'.format(unit=angUnit),fontsize=14)
		plt.xlabel(r'$\theta$ [{unit}]'.format(unit=angUnit),fontsize=14)
		plt.title('Theta vs Phi',fontsize=16)
		pltName = 'plot_theta_vs_phi.png'
		plt.savefig(pltName, dpi=300, bbox_inches='tight')
		print(f'Created {pltName}')
		

	if procedure == 'thvsb':
		# Place values on variables
		theta_in = column_lists['theta']
		b_out = column_lists['b_out']

		# b vs theta
	
		plt.figure(figsize=(8,6), facecolor='w')
		plt.plot(theta_in, b_out)
		plt.xlim(min(theta_in), max(theta_in))
		plt.ylabel(r'$b$ [fm]',fontsize=14)
		plt.xlabel(r'$\theta$ [{unit}]'.format(unit=angUnit),fontsize=14)
		plt.title('Impact parameter as function of the scattering angle',fontsize=16)
		pltName = 'plot_b_vs_theta.png'
		plt.savefig(pltName, dpi=300, bbox_inches='tight')
		print(f'Created {pltName}')
	

		# theta vs b

		plt.figure(figsize=(8,6), facecolor='w')
		plt.plot(b_out, theta_in)
		plt.xlim(min(b_out), max(b_out))
		plt.ylabel(r'$\theta$ [{unit}]'.format(unit=angUnit),fontsize=14)
		plt.xlabel(r'$b$ [fm]',fontsize=14)
		plt.title('Scattering angle as function of the impact parameter',fontsize=16)
		pltName = 'plot_theta_vs_b.png'
		plt.savefig(pltName, dpi=300, bbox_inches='tight')
		print(f'Created {pltName}')

	if procedure == 'bvsd':
		b_in = column_lists['b_in']
		closestDistance = column_lists['closestDistance']
		
		plt.figure(figsize=(8,6), facecolor='w')
		plt.plot(b_in, closestDistance)
		plt.xlim(min(b_in), max(b_in))
		plt.ylabel(r'Minimum distance'.format(unit=angUnit),fontsize=14)
		plt.xlabel(r'$b$',fontsize=14)
		plt.title('Closest distance as function of the impact parameter',fontsize=16)
		pltName = 'plot_closest_vs_b.png'
		plt.savefig(pltName, dpi=300, bbox_inches='tight')
		print(f'Created {pltName}')


	# Differential cross section 
	if procedure == 'xsec':
	
		theta_in = column_lists['theta']
		
		# Place values on dictionary for each DifCrossSec correction
		difCrossSec = {}
		difCrossSec['Rutherford'] = column_lists['difCrossSec_Ruth']
		if mott:
			difCrossSec['Mott'] = column_lists['difCrossSec_Mott']
		if recoil:
			difCrossSec['Target Recoil'] = column_lists['difCrossSec_Recoil']
		if diracProton:
			difCrossSec['Dirac Proton'] = column_lists['difCrossSec_diracProton']
		if formFactor:
			difCrossSec['Form Factor'] = column_lists['difCrossSec_formFactor']
		if rosenbluth:
			difCrossSec['Rosenbluth'] = column_lists['difCrossSec_rosenbluth']

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
			
			# Extract the x values from the JSON data
			hofAngles = [float(entry["x"][0]["value"]) for entry in real["values"]]

			# Get energy value for labels
			col_names = [E["value"] for E in real["qualifiers"]['E']]


		# Plotting of Hoftstadter data
		if hof25:
			if kinEn != 25e6:
				print(f"Chosen kinetic energy {int(kinEn*1e-6)}MeV is not 25 MeV, therefore data was not plotted.")
			else:
				# Gets values from file
				hof_difCrossSec_25 = []
				hofAngles25 = []
				with open(sys.path[0] + '/../data/hoftstadter25.csv', 'r') as hof25:
					for line in hof25:
						angle, value = line.strip().split(';')
						hofAngles25.append(float(angle))
						hof_difCrossSec_25.append(float(value))

				# Converts data units
				hofAngles25, hof_difCrossSec_25= convert_data(hofAngles25, hof_difCrossSec_25, None, cross_section_variable)

				# Normalizes simluation
				for key in difCrossSec:
					difCrossSec[key] = normalize_sim(hofAngles25, hof_difCrossSec_25, theta_in, difCrossSec[key])

				plt.scatter(hofAngles25, hof_difCrossSec_25, label="Hoftstadter 25 Mev") 


		if hof125:
			if kinEn != 125e6:
				print(f"Chosen kinetic energy {int(kinEn*1e-6)}MeV is not 125 MeV, therefore data was not plotted.")
			else:
				# Gets values from file
				hof_difCrossSec_125 = []
				hofAngles125 = []
				with open(sys.path[0] + '/../data/hoftstadter125.csv', 'r') as hof125:
					for line in hof125:
						angle, value = line.strip().split(';')
						hofAngles125.append(float(angle))
						hof_difCrossSec_125.append(float(value))
			
				# Converts data units
				hofAngles125, hof_difCrossSec_125= convert_data(hofAngles125, hof_difCrossSec_125, None, cross_section_variable)
			
				# Normalizes simluation
				for key in difCrossSec:
					difCrossSec[key] = normalize_sim(hofAngles125, hof_difCrossSec_125, theta_in, difCrossSec[key])
				
				plt.scatter(hofAngles125, hof_difCrossSec_125, marker='s', color='black',  label="Hoftstadter 125 Mev") 


		if hof300:
			if kinEn != 300e6:
				print(f"Chosen kinetic energy {int(kinEn*1e-6)}MeV is not 300 MeV, therefore data was not plotted.")
			else:
				# Gets values from json
				hof_difCrossSec_300 = [float(entry["y"][0]["value"]) for entry in real["values"]]
				yerr = [float(entry["y"][2]["errors"][0]["symerror"]) for entry in real["values"]]

				# Converts data units
				hofAngles, hof_difCrossSec_300, yerr = convert_data(hofAngles, hof_difCrossSec_300, yerr, cross_section_variable)

				# Normalizes simluation
				for key in difCrossSec:
					difCrossSec[key] = normalize_sim(hofAngles, hof_difCrossSec_300, theta_in, difCrossSec[key])

				plt.errorbar(hofAngles, hof_difCrossSec_300, yerr, capsize = 3, ls='none', label="Hoftstadter "+col_names[0]) 
		

		if hof400:
			if kinEn != 400e6:
				print(f"Chosen kinetic energy {int(kinEn*1e-6)}MeV is not 400 MeV, therefore data was not plotted.")
			else:
				# Gets values from json
				hof_difCrossSec_400 = [float(entry["y"][1]["value"]) for entry in real["values"]]
				yerr = [float(entry["y"][2]["errors"][0]["symerror"]) for entry in real["values"]]

				# Converts data units
				hofAngles, hof_difCrossSec_400, yerr = convert_data(hofAngles, hof_difCrossSec_400, yerr, cross_section_variable)

				# Normalizes simluation
				for key in difCrossSec:
					difCrossSec[key] = normalize_sim(hofAngles, hof_difCrossSec_400, theta_in, difCrossSec[key])

				plt.errorbar(hofAngles, hof_difCrossSec_400, yerr, capsize = 3, ls='none', label="Hoftstadter "+col_names[1])
	

		if hof550:
			if kinEn != 550e6:
				print(f"Chosen kinetic energy {int(kinEn*1e-6)}MeV is not 550 MeV!, therefore data was not plotted.")
			else:
				# Gets values from json
				hof_difCrossSec_550 = [float(entry["y"][3]["value"]) for entry in real["values"]]
				yerr = [float(entry["y"][3]["errors"][0]["symerror"]) for entry in real["values"]]

				# Converts data units
				hofAngles, hof_difCrossSec_550, yerr = convert_data(hofAngles, hof_difCrossSec_550, yerr, cross_section_variable)

				# Normalizes simluation
				for key in difCrossSec:
					difCrossSec[key] = normalize_sim(hofAngles, hof_difCrossSec_550, theta_in, difCrossSec[key])

				plt.errorbar(hofAngles, hof_difCrossSec_550, yerr, capsize = 3, ls='none', label="Hoftstadter "+col_names[3])


		if any([hof25, hof125, hof300, hof400, hof550]):
			pltName += '_hoftstadter'


		if geiger:
			if kinEn != 125e6:
				print(f"Chosen kinetic energy {int(kinEn*1e-6)}MeV is not 125 MeV!, therefore data was not plotted.")
			else:
				# Gets values from file
				geiger_difCrossSec_125 = []
				geigerAngles125 = []
				with open(sys.path[0] + '/../data/GeigerMarsden.csv', 'r') as hof125:
					for line in hof125:
						angle, value = line.strip().split(';')
						geigerAngles125.append(float(angle))
						geiger_difCrossSec_125.append(float(value))
	
				# Converts data units
				geigerAngles125, geiger_difCrossSec_125 = convert_data(geigerAngles125, geiger_difCrossSec_125, None, cross_section_variable)
			
				# Normalizes simluation
				for key in difCrossSec:
					difCrossSec[key] = normalize_sim(geigerAngles125, geiger_difCrossSec_125, theta_in, difCrossSec[key])
				
				plt.scatter(geigerAngles125, geiger_difCrossSec_125, color='black',  label="GeigerMarsden") 

				pltName += '_GeigerMarsden'


		# Plots every Differential cross section
		for key in difCrossSec:
			plt.plot(theta_in, difCrossSec[key], label=key)
		
		# Finalizes plot
		plt.xlim(min(theta_in), max(theta_in))
		pltName += f'_{int(kinEn*1e-6)}MeV'
		plt.legend()
		plt.savefig(pltName, dpi=300, bbox_inches='tight')
		print(f'Created {pltName}.png')


if __name__ == '__main__':
    plot('output.dat','input.dat')