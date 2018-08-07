from metalicity import *

if __name__ == "__main__":
	tag = 1
	sca = 0
	#nline = 42
	#ncore = 6
	rep0 = 'halo1_jj/'
	ldir = ['NL4_zoom_wdm/'+rep0, 'NL4_zoom_cdm/'+rep0]
	#ldir = ['halo1_jj_wdm/','halo1_jj_cdm/']
	Zsh = 1e-4

	out0 = []
	out1 = []
	if tag==0:
		for sn in range(20,21):
			out0.append(metal(sn=sn, indm=1, rep=ldir[0]))
			out1.append(metal(sn=sn, indm=0, rep=ldir[1]))
		data = np.hstack([out0,out1])
		data = data.T
		totxt(rep0+'metalicity.txt',data,0,0,0)
	else:
		data = retxt(rep0+'metalicity.txt',6,0,0)

	plt.figure()
	plt.plot(data[2], data[0], label=lmodel[1],marker='^')
	plt.plot(data[2], data[3], label=lmodel[0],ls='--',marker='o')
	plt.xlabel(r'$z$')
	plt.ylabel(r'$\bar{Z}\ [Z_{\odot}]$')
	plt.legend()
	if sca==0:
		plt.yscale('log')
	plt.tight_layout()
	if sca==0:
		plt.savefig(rep0+'logZ_z.pdf')
	else:
		plt.savefig(rep0+'Zbar_z.pdf')	

	plt.figure()
	plt.plot(data[2], data[1], label=lmodel[1],marker='^')
	plt.plot(data[2], data[4], label=lmodel[0],ls='--',marker='o')
	plt.xlabel(r'$z$')
	plt.ylabel(r'$f_{\mathrm{PopII}}$')
	plt.legend()
	if sca==0:
		plt.yscale('log')
	plt.tight_layout()
	if sca==0:
		plt.savefig(rep0+'logfrac_z.pdf')
	else:
		plt.savefig(rep0+'frac_z.pdf')
