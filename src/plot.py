from card_reader import read_card, _raise_missing_card_error 
import matplotlib.pyplot as plt
import numpy as np
import json
import sys

def plot(output,card):

	# Read file
	data = np.loadtxt(output, delimiter=",")

	#Read parameter input card 
	parameters = read_card(card)
	angUnit = parameters[3]
	mott = parameters[6]
	
	theta_in = data[:, 0]  
	b_out = data[:, 1]
	dsig_dtheta = data[:, 2]
	if mott == "true":
		dsig_dtheta_Mott = data[:, 3]
	

	# Either cos, theta or omega. 
	var = 'theta'

	# Set Hofstadter True or False
	hof = True

	# Read Hofstadter data 
	with open(sys.path[0] + '/../data/Hofstadter.json', 'r') as json_file:
	    real = json.load(json_file)


	# b vs theta

	plt.figure(figsize=(8,6), facecolor='w')
	plt.plot(theta_in, b_out)
	plt.ylabel(r'$b$ [fm]',fontsize=14)
	plt.xlabel(r'$\theta$ [{unit}]'.format(unit=angUnit),fontsize=14)
	plt.title('Impact parameter as function of the scattering angle',fontsize=16)

	plt.savefig('plot_b_vs_theta.png', dpi=300, bbox_inches='tight')


	# theta vs b

	plt.figure(figsize=(8,6), facecolor='w')
	plt.plot(b_out, theta_in)
	plt.ylabel(r'$\theta$ [{unit}]'.format(unit=angUnit),fontsize=14)
	plt.xlabel(r'$b$ [fm]',fontsize=14)
	plt.title('Scattering angle as function of the impact parameter',fontsize=16)

	plt.savefig('plot_theta_vs_b.png', dpi=300, bbox_inches='tight')


	# differential cross section vs dtheta/dcos(theta)/domega

	pltName='plot_dsig_'
	plt.figure(figsize=(8,6), facecolor='w')
	plt.plot(theta_in, dsig_dtheta, label='Rutherford')
	if mott == "true":
		plt.plot(theta_in, dsig_dtheta_Mott, label='Mott')
	plt.yscale("log")

	# Change labels
	if var == 'cos':
		pltName += 'dcostheta_vs_costheta'
		plt.xlabel(r'$cos(\theta)$'.format(unit=angUnit),fontsize=14)
		plt.ylabel(r'$d\sigma/dcos(\theta)$',fontsize=14)
		plt.title(r'Distribution of $d\sigma/dcos(\theta)$ as function of the scattering angle',fontsize=16)

	elif var == 'theta': 
		pltName += 'dtheta_vs_theta'
		plt.xlabel(r'$\theta$ [{unit}]'.format(unit=angUnit),fontsize=14)
		plt.ylabel(r'$d\sigma/d\theta$',fontsize=14)
		plt.title(r'Distribution of $d\sigma/d\theta$ as function of the scattering angle',fontsize=16)
		#Use hoftstadter data
		if hof:
			pltName += '_hoftstadter'
			# Extract the x and y values from the JSON data
			hofAngles = [float(entry["x"][0]["value"]) for entry in real["values"]]
			dsigdOmega_300 = [float(entry["y"][0]["value"]) for entry in real["values"]]
			dsigdOmega_400 = [float(entry["y"][1]["value"]) for entry in real["values"]]
			dsigdOmega_500 = [float(entry["y"][2]["value"]) if entry["y"][2]["value"] != "-" else None for entry in real["values"]]
			dsigdOmega_550 = [float(entry["y"][3]["value"]) for entry in real["values"]]

			plt.plot(hofAngles, dsigdOmega_300, label='Hofstadter 300 MeV')
			plt.plot(hofAngles, dsigdOmega_400, label='Hofstadter 400 MeV')
			plt.plot(hofAngles, dsigdOmega_500, label='Hofstadter 500 MeV')
			plt.plot(hofAngles, dsigdOmega_550, label='Hofstadter 550 MeV')
		
	elif var == 'omega':
		pltName += 'domega_vs_theta'
		plt.xlabel(r'$\theta$ [{unit}]'.format(unit=angUnit),fontsize=14)
		plt.ylabel(r'$d\sigma/d\Omega$',fontsize=14)
		plt.title(r'Distribution of $d\sigma/d\Omega$ as function of the scattering angle',fontsize=16)
	plt.legend()
	plt.savefig(pltName+',png', dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    plot('output.dat','input.dat')
