from card_reader import read_card, _raise_missing_card_error 
import matplotlib.pyplot as plt
import numpy as np
import json
import sys

def plot(output,card_name):

	# Read file
	data = np.loadtxt(output, delimiter=",")

	#Read parameter input card 
	parameters = read_card(card_name)
	kinEn = parameters[0]
	angUnit = parameters[3]
	mott = parameters[6]
	hof300 = parameters[8]
	hof400 = parameters[9]
	hof500 = parameters[10]
	hof550 = parameters[11]
	
	theta_in = data[:, 0]  
	b_out = data[:, 1]
	dsig_dtheta = data[:, 2]
	if mott == "true":
		dsig_dtheta_Mott = data[:, 3]
	

	# Either cos, theta or omega. 
	var = 'theta'


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
		"""
		Hoftstadter data
		Data from 500MeV has a datapoint with zero uncertainty, so it is ignored
		"""
		# Extract the x and y values from the JSON data
		hofAngles = [float(entry["x"][0]["value"]) for entry in real["values"]]
		# Get energy value for labels
		col_names = [E["value"] for E in real["qualifiers"]['E']]

		if hof300 == 'true':
			hof = True
			dsigdOmega_300 = [float(entry["y"][0]["value"]) for entry in real["values"]]
			plt.errorbar(hofAngles, dsigdOmega_300, yerr=[float(entry["y"][0]["errors"][0]["symerror"]) for entry in real["values"]], capsize = 3, ls='none', label="Hoftstadter "+col_names[0]) 
			
		if hof400 == 'true':
			hof = True
			dsigdOmega_400 = [float(entry["y"][1]["value"]) for entry in real["values"]]
			plt.errorbar(hofAngles, dsigdOmega_400, yerr=[float(entry["y"][1]["errors"][0]["symerror"]) for entry in real["values"]],  capsize = 3, ls='none', label="Hoftstadter "+col_names[1])

		if hof500 == 'true':
			hof = True
			dsigdOmega_500 = [float(entry["y"][2]["value"]) for entry in real["values"] if entry["y"][2]["value"] != "-"]
			plt.errorbar(hofAngles[1:], dsigdOmega_500, yerr=[float(entry["y"][2]["errors"][0]["symerror"]) for entry in real["values"] if entry["y"][2]["errors"][0]["symerror"] != 0], capsize = 3, ls='none', label="Hoftstadter "+col_names[2])

		if hof550 == 'true':
			hof = True
			dsigdOmega_550 = [float(entry["y"][3]["value"]) for entry in real["values"]]
			plt.errorbar(hofAngles, dsigdOmega_550, yerr=[float(entry["y"][3]["errors"][0]["symerror"]) for entry in real["values"]],  capsize = 3, ls='none', label="Hoftstadter "+col_names[3])
		
		if hof == True:
			pltName += '_hoftstadter'
	elif var == 'omega':
		pltName += 'domega_vs_theta'
		plt.xlabel(r'$\theta$ [{unit}]'.format(unit=angUnit),fontsize=14)
		plt.ylabel(r'$d\sigma/d\Omega$',fontsize=14)
		plt.title(r'Distribution of $d\sigma/d\Omega$ as function of the scattering angle',fontsize=16)
	plt.legend()
	plt.savefig(pltName, dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    plot('output.dat','input.dat')