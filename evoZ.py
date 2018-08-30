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

	sn0 = 25
	sn1 = 25

	out0 = []
	out1 = []
	if tag==0:
		for sn in range(0,26):
			if sn<=sn0:
				out0.append(metal(sn=sn, indm=1, rep=ldir[0]))
			if sn<=sn1:
				out1.append(metal(sn=sn, indm=0, rep=ldir[1]))
		#data = np.hstack([out0,out1])
		#data = data.T
		out0 = np.array(out0).T
		out1 = np.array(out1).T
		totxt(rep0+'metalicity'+'_'+lmodel[1]+'.txt',out0,0,0,0)
		totxt(rep0+'metalicity'+'_'+lmodel[0]+'.txt',out1,0,0,0)
	else:
		out0 = np.array(retxt(rep0+'metalicity'+'_'+lmodel[1]+'.txt', 6, 0, 0))
		out1 = np.array(retxt(rep0+'metalicity'+'_'+lmodel[0]+'.txt', 6, 0, 0))

	plt.figure()
	plt.plot(out0[2][out0[0]>0], out0[0][out0[0]>0], label='Total, '+lmodel[1],lw=2)#,marker='^')
	plt.plot(out1[2][out1[0]>0], out1[0][out1[0]>0], label='Total, '+lmodel[0],ls='--',lw=2)#,marker='^')
	plt.plot(out0[2][out0[4]>0], out0[4][out0[4]>0], label='PopII, '+lmodel[1],marker='^',lw=1)
	plt.plot(out1[2][out1[4]>0], out1[4][out1[4]>0], label='PopII, '+lmodel[0],ls='--',marker='^',lw=1)
	plt.plot(out0[2][out0[5]>0], out0[5][out0[5]>0], label='PopIII, '+lmodel[1],marker='o',lw=1)
	plt.plot(out1[2][out1[5]>0], out1[5][out1[5]>0], label='PopIII, '+lmodel[0],ls='--',marker='o',lw=1)
	plt.xlabel(r'$z$')
	plt.ylabel(r'$\langle Z\rangle\ [Z_{\odot}]$')
	plt.legend()
	if sca==0:
		plt.yscale('log')
	plt.tight_layout()
	if sca==0:
		plt.savefig(rep0+'logZ_z.pdf')
	else:
		plt.savefig(rep0+'Zbar_z.pdf')	

	plt.figure()
	plt.plot(out0[2][out0[1]>0], out0[1][out0[1]>0], label='Volume, '+lmodel[1],marker='^')
	plt.plot(out1[2][out1[1]>0], out1[1][out1[1]>0], label='Volume, '+lmodel[0],ls='--',marker='^')
	plt.plot(out0[2][out0[3]>0], out0[3][out0[3]>0], label='Mass, '+lmodel[1],marker='o')
	plt.plot(out1[2][out1[3]>0], out1[3][out1[3]>0], label='Mass, '+lmodel[0],ls='--',marker='o')
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
