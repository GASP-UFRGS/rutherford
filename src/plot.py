import matplotlib.pyplot as plt
import numpy as np

def plot(output):

	# Read file
	data = np.loadtxt(output, delimiter=",")

	# Get angle unit
	with open(output) as file:
	    lines = file.readlines()
	    angUnit = lines[-1].lstrip("# ").strip()

	theta_in = data[:, 0]  
	b_out = data[:, 1]
	dsig_dtheta = data[:, 2]
	dsig_dtheta_Mott = data[:, 3]
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
	plt.plot(theta_in, dsig_dtheta_Mott, label='Mott')
	plt.yscale("log")
	plt.xlabel(r'$\theta$ [{unit}]'.format(unit=angUnit),fontsize=14)
	plt.ylabel(r'$d\sigma/dcos(\theta)$',fontsize=14)
	plt.title(r'Distribution of $d\sigma/dcos(\theta)$ as function of the scattering angle',fontsize=16)
	plt.legend()

	plt.savefig('plot_dsig_dtheta_vs_theta.png', dpi=300, bbox_inches='tight')



if __name__ == '__main__':
    plot('output.dat')
