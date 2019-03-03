from radio import *
Vz = 4**3 * 0.3783187*0.3497985*0.399293243 /0.6774**3

if __name__ == "__main__":
	sca =int(sys.argv[1])
	rep0 = 'halo1_jj/'
	#rep0 = 'halo1_jj_new/'
	ldir = ['NL4_zoom_wdm/'+rep0, 'NL4_zoom_cdm/'+rep0]
	#ldir = ['halo1_jj_wdm/', 'halo1_jj_cdm/']
	#ldir = ['halo1_jj_wdm_new/','halo1_jj_cdm_new/']
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

	print('Volume of the zoom-in region: {} [Mpc^3]'.format(Vz))
	print('SFR_popIII ratio: {}, SFR_popII ratio: {}'.format(10**fs0(8.5)/10**fs1(8.5), 10**fs0_(8.5)/10**fs1_(8.5)))

	repref = 'boyuan_sfrd/'
	refIII = np.array(retxt(repref+'popIII_sfr_s2_skx.txt',4,0,0))
	refII = np.array(retxt(repref+'popII_sfr_s2_skx.txt',4,0,0))

	down1, up1 = 1e-6, 10
	plt.figure()
	#plt.plot(lz0[d0_III[2]>0],(d0_III[2]+d0_II[2])[d0_III[2]>0],label=lmodel[1])
	plt.plot(lz0[d0_III[2]>0],d0_III[2][d0_III[2]>0]/Vz,label='PopIII, '+lmodel_[1],lw=3)
	plt.plot(lz1[d1_III[2]>0],d1_III[2][d1_III[2]>0]/Vz,label='PopIII, '+lmodel_[0])
	plt.plot(lz0_[d0_II[2]>0],d0_II[2][d0_II[2]>0]/Vz,label='PopII, '+lmodel_[1],ls='--',lw=3)
	plt.plot(lz1_[d1_II[2]>0],d1_II[2][d1_II[2]>0]/Vz,label='PopII, '+lmodel_[0],ls='--')
	plt.plot(1/refIII[0][refIII[3]>0]-1, refIII[3][refIII[3]>0],label='PopIII, Jaacks et al. (2018a)',lw=0.5)
	plt.plot(1/refII[0][refII[3]>0]-1, refII[3][refII[3]>0],label='PopII, Jaacks et al. (2018a)',ls='--',lw=0.5)
	#plt.plot(lz1[d1_III[2]>0],(d1_III[2]+d1_II[2])[d1_III[2]>0],label=lmodel_[0],ls='--')
	plt.fill_between([16, 19],[up1, up1],[down1, down1],label='EDGES',facecolor='gray', alpha=0.5)
	plt.text(8,2e-6,r'$V_{Z,c}\simeq 11\ \mathrm{Mpc^{3}}$')
	plt.ylim(down1, up1)
	plt.xlim(6, 27)
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

	fs0 = interp1d(np.log10([TZ(z)/1e6/YR for z in lz0]), np.log10(d0_III[2]+Vz*1e-8*(d0_III[2]==0)))
	fs0_ = interp1d(np.log10([TZ(z)/1e6/YR for z in lz0]), np.log10(d0_II[2]+Vz*1e-8*(d0_II[2]==0)))
	fs1 = interp1d(np.log10([TZ(z)/1e6/YR for z in lz1]), np.log10(d1_III[2]+Vz*1e-8*(d1_III[2]==0)))
	fs1_ = interp1d(np.log10([TZ(z)/1e6/YR for z in lz1]), np.log10(d1_II[2]+Vz*1e-8*(d1_II[2]==0)))

	def cumsfr(f, lz, tau = 3, n = 10):
		lt = [TZ(z)/1e6/YR for z in lz]
		t0 = np.min(lt)
		lt0 = [max(t-tau, t0) for t in lt]
		out = np.zeros(len(lt))
		for i in range(len(lt)):
			dt = np.linspace(lt0[i], lt[i], n)
			out[i] = np.trapz(10**f(np.log10(dt)), dt)*1e6
		return out
		
	lsm0 = cumsfr(fs0, lz0, 3)  #np.cumsum(d0_III[1]*d0_III[2])
	lsm0_ = cumsfr(fs0_, lz0, 10) #np.cumsum(d0_II[1]*d0_II[2])
	lsm1 = cumsfr(fs1, lz1, 3) #np.cumsum(d1_III[1]*d1_III[2])
	lsm1_ = cumsfr(fs1_, lz1, 10) #np.cumsum(d1_II[1]*d1_II[2])
	for i in range(len(lsm1_)):
		if lsm1_[i]==0 and lz1_[i]<23:
			lsm1_[i] = 1e-4

	f0 = interp1d(lz0[lsm0>0],np.log10(lsm0[lsm0>0]))
	f0_ = interp1d(lz0_[lsm0_>0],np.log10(lsm0_[lsm0_>0]))
	f1 = interp1d(lz1[lsm1>0],np.log10(lsm1[lsm1>0]))
	f1_ = interp1d(lz1_[lsm1_>0],np.log10(lsm1_[lsm1_>0]))

	down2, up2 = 5e1/Vz,2e7/Vz
	plt.figure()
	plt.plot(lz0[lsm0>0],lsm0[lsm0>0]/Vz,label='PopIII, '+lmodel_[1],lw=3)
	plt.plot(lz1[lsm1>0],lsm1[lsm1>0]/Vz,label='PopIII, '+lmodel_[0])
	plt.plot(lz0_[lsm0_>0],lsm0_[lsm0_>0]/Vz,label='PopII, '+lmodel_[1],ls='--',lw=3)
	plt.plot(lz1_[lsm1_>0],lsm1_[lsm1_>0]/Vz,label='PopII, '+lmodel_[0],ls='--')
	plt.fill_between([16, 19],[up2, up2],[down2, down2],label='EDGES',facecolor='gray', alpha=0.5)
	plt.text(6.5,1e1,r'$V_{Z,c}\simeq 11\ \mathrm{Mpc^{3}}$')
	plt.xlabel(r'$z$')
	plt.ylabel(r'$\rho_{\star}\ [M_{\odot}\ \mathrm{Mpc^{-3}}]$')
	plt.ylim(down2, up2)
	plt.xlim(6, 27)
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
