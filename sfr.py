from radio import *
Vz = 4**3 * 0.3783187*0.3497985*0.399293243 /0.6774**3

if __name__ == "__main__":
	sca =int(sys.argv[1])
	rep0 = 'halo1_jj/'
	#ldir = ['NL4_zoom_wdm/'+rep0, 'NL4_zoom_cdm/'+rep0]
	ldir = ['halo1_jj_wdm/', 'halo1_jj_cdm/']
	d0_III = np.array(retxt(ldir[0]+'popIII_sfr.txt',4,0,0))
	d0_II = np.array(retxt(ldir[0]+'popII_sfr.txt',4,0,0))
	d1_III = np.array(retxt(ldir[1]+'popIII_sfr.txt',4,0,0))
	d1_II = np.array(retxt(ldir[1]+'popII_sfr.txt',4,0,0))
	lz0 = (1/d0_III[0]-1)
	lz0_ = (1/d0_II[0]-1)
	lz1 = (1/d1_III[0]-1)
	lz1_ = (1/d1_II[0]-1)

	fs0 = interp1d(lz0[d0_III[2]>0],np.log10(d0_III[2][d0_III[2]>0]))
	fs0_ = interp1d(lz0_[d0_II[2]>0],np.log10(d0_II[2][d0_II[2]>0]))
	fs1 = interp1d(lz1[d1_III[2]>0],np.log10(d1_III[2][d1_III[2]>0]))
	fs1_ = interp1d(lz1_[d1_II[2]>0],np.log10(d1_II[2][d1_II[2]>0]))

	print('SFR_popIII ratio: {}, SFR_popII ratio: {}'.format(10**fs0(8.5)/10**fs1(8.5), 10**fs0_(8.5)/10**fs1_(8.5)))

	plt.figure()
	#plt.plot(lz0[d0_III[2]>0],(d0_III[2]+d0_II[2])[d0_III[2]>0],label=lmodel[1])
	plt.plot(lz0[d0_III[2]>0],d0_III[2][d0_III[2]>0]/Vz,label='PopIII, '+lmodel[1])
	plt.plot(lz1[d1_III[2]>0],d1_III[2][d1_III[2]>0]/Vz,label='PopIII, '+lmodel[0],ls='--')
	plt.plot(lz0_[d0_II[2]>0],d0_II[2][d0_II[2]>0]/Vz,label='PopII, '+lmodel[1],ls='-.')
	plt.plot(lz1_[d1_II[2]>0],d1_II[2][d1_II[2]>0]/Vz,label='PopII, '+lmodel[0],ls=':')
	#plt.plot(lz1[d1_III[2]>0],(d1_III[2]+d1_II[2])[d1_III[2]>0],label=lmodel[0],ls='--')
	plt.xlabel(r'$z$')
	plt.ylabel(r'$\mathrm{SFRD}\ [M_{\odot}\ \mathrm{yr^{-1}\ Mpc^{-3}}]$')
	if sca>0:
		plt.yscale('log')
	plt.legend()
	plt.tight_layout()
	if sca==0:
		plt.savefig(rep0+'SFRD_z.pdf')
	else:
		plt.savefig(rep0+'logSFRD_z.pdf')

	lsm0 = np.cumsum(d0_III[1]*d0_III[2])
	lsm0_ = np.cumsum(d0_II[1]*d0_II[2])
	lsm1 = np.cumsum(d1_III[1]*d1_III[2])
	lsm1_ = np.cumsum(d1_II[1]*d1_II[2])

	f0 = interp1d(lz0[lsm0>0],np.log10(lsm0[lsm0>0]))
	f0_ = interp1d(lz0_[lsm0_>0],np.log10(lsm0_[lsm0_>0]))
	f1 = interp1d(lz1[lsm1>0],np.log10(lsm1[lsm1>0]))
	f1_ = interp1d(lz1_[lsm1_>0],np.log10(lsm1_[lsm1_>0]))

	plt.figure()
	plt.plot(lz0[lsm0>0],lsm0[lsm0>0],label='PopIII, '+lmodel[1])
	plt.plot(lz1[lsm1>0],lsm1[lsm1>0],label='PopIII, '+lmodel[0],ls='--')
	plt.plot(lz0_[lsm0_>0],lsm0_[lsm0_>0],label='PopII, '+lmodel[1],ls='-.')
	plt.plot(lz1_[lsm1_>0],lsm1_[lsm1_>0],label='PopII, '+lmodel[0],ls=':')
	plt.xlabel(r'$z$')
	plt.ylabel(r'$M_{*}\ [M_{\odot}]$')
	if sca>0:
		plt.yscale('log')
	plt.legend()
	plt.tight_layout()
	if sca==0:
		plt.savefig(rep0+'Mstar_z.pdf')
	else:
		plt.savefig(rep0+'logMstar_z.pdf')
	#plt.show()

	print('M_popIII ratio: {}, M_popII ratio: {}'.format(10**f0(8.5)/10**f1(8.5), 10**f0_(8.5)/10**f1_(8.5)))
