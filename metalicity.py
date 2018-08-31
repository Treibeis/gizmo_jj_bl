from radio import *

def metal(sn = 50, rep = './', indm = 0, edge = [0.01, 100.0], base = 'snapshot', ext = '.hdf5', Zth = 1e-4, mode=0, Zsun = 0.0134):
	ds = yt.load(rep+base+'_'+str(sn).zfill(3)+ext)
	z = ds['Redshift']
	ad = ds.all_data()
	keys = ds.field_list
	tag = np.sum([x[0] == 'PartType4' for x in keys])
	Ngas = len(ad[('PartType0', 'Metallicity_00')])
	Zraw = ad[('PartType0', 'Metallicity_00')]/Zsun#+ad[('PartType0', 'Metallicity_01')]+ad[('PartType0', 'Metallicity_02')])
	Zraw0 = ad[('PartType0', 'Metallicity_01')]/Zsun
	Zraw1 = ad[('PartType0', 'Metallicity_02')]/Zsun
	lV = ad[('PartType0', 'Masses')]/ad[('PartType0', 'Density')]
	popII = Zraw > Zth
	lZ = Zraw[popII]
	#if tag>0:
	#	Nstar = len(ad[('PartType4', 'Metallicity_00')])
	#	Zraw_ = (ad[('PartType4', 'Metallicity_00')]+ad[('PartType4', 'Metallicity_01')]+ad[('PartType4', 'Metallicity_02')])
	#	popII_ = Zraw_ > Zth
	#	lZ_ = Zraw_[popII_]
	
	if mode>0:
		nd = ad[('PartType0','density')].to_equivalent("cm**-3", "number_density",mu=mmw(ad[('PartType0','Primordial HII')]))
		rn = np.percentile(np.log10(nd), edge)
		#rZ = np.percentile(np.log10(Zraw), edge)
		#rZ[0] = max(rZ[0], -5)
		rZ = [-5, 0]
		plt.figure()
		plt.subplot(111)
		plt.hist2d(np.log10(nd),np.log10(Zraw),bins=100,norm=LogNorm(),range=[rn,rZ])
		cb = plt.colorbar()
		cb.set_label(r'$\log(N)$')
		cb.set_clim(1.0,6e2)
		plt.plot(rn,[-4,-4],'k--',label=r'$Z_{\mathrm{crit}}$')
		plt.legend()
		plt.ylabel(r'$\log(Z\ [Z_{\odot}])$')
		plt.xlabel(r'$\log(n\ [\mathrm{cm^{-3}}])$')
		#plt.title(r'$[\mathrm{H_{2}/H}]-T$ phase diagram for '+lmodel[indm]+' at $z=$'+str(int(ds['Redshift']*100)/100),size=12)
		plt.tight_layout()
		plt.savefig(rep+'Zn_'+lmodel[indm]+'_'+str(sn)+'.pdf')

	#if tag>0:
	#	rat = (len(lZ)+len(lZ_))/(Ngas+Nstar)
	#else:
	rat0 = np.array(np.sum(lV[popII])/np.sum(lV))
	rat1 = np.array(np.sum(ad[('PartType0', 'Masses')][popII])/np.sum(ad[('PartType0', 'Masses')]))
	Zbar = np.average(Zraw)
	Zbar0 = np.average(Zraw0)
	Zbar1 = np.average(Zraw1)
	return np.array([Zbar, rat0, z, rat1, Zbar0, Zbar1])

if __name__ == "__main__":
	indm = int(sys.argv[2])
	sn = int(sys.argv[1])
	if indm==0:
		rep0 = 'halo1_jj_cdm/'
	else:
		rep0 = 'halo1_jj_wdm/'
	out = metal(sn=sn, indm=indm, rep=rep0,mode=1)
	print('Average Z: {}, Mass fraction: {}'.format(out[0], out[1]))
