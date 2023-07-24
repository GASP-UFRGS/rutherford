from card_reader import read_card, _raise_missing_card_error 
import matplotlib.pyplot as plt
import numpy as np

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
	

	# Triggers for changing label to dcos(theta)
	cos = False
	
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

	# theta vs dsig_dtheta

	plt.figure(figsize=(8,6), facecolor='w')
	plt.plot(theta_in, dsig_dtheta, label='Rutherford')
	if mott == "true":
		plt.plot(theta_in, dsig_dtheta_Mott, label='Mott')
	plt.yscale("log")
	plt.legend()
	if cos:
		plt.xlabel(r'$cos(\theta)$ [{unit}]'.format(unit=angUnit),fontsize=14)
		plt.ylabel(r'$d\sigma/dcos(\theta)$',fontsize=14)
		plt.title(r'Distribution of $d\sigma/dcos(\theta)$ as function of the scattering angle',fontsize=16)
		plt.savefig('plot_dsig_dtheta_vs_costheta.png', dpi=300, bbox_inches='tight')
	else: 
		plt.xlabel(r'$\theta$ [{unit}]'.format(unit=angUnit),fontsize=14)
		plt.ylabel(r'$d\sigma/d\theta$',fontsize=14)
		plt.title(r'Distribution of $d\sigma/d\theta$ as function of the scattering angle',fontsize=16)
		plt.savefig('plot_dsig_dtheta_vs_theta.png', dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    plot('output.dat','input.dat')
