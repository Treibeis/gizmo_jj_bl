
from radio import *

if __name__ == "__main__":
	ncore = 4
	Tsh = 1e4
	nsh = 1e8
	indm = 1
	sn = int(sys.argv[1])
	mode = int(sys.argv[2])
	bins = int(sys.argv[3])
	#nn = int(sys.argv[3])
	if len(sys.argv)>5:
		low, up = int(sys.argv[4]), int(sys.argv[5])
	else:
		low, up = 1900, 2000
	if mode==0:
		mesh = mesh_3d(low, up, bins, low, up, bins, low, up, bins)
		d0 = grid(mesh,sn,box =[[low]*3,[up]*3],ncore=ncore,Tsh=Tsh,nsh=nsh)
		#d1 = grid(mesh,sn,box =[[low]*3,[up]*3],mode=1,ncore=ncore,Tsh=Tsh)
	else:
		d0 = read_grid(sn,nb = bins)
		#d1 = read_grid(sn,nb = bins, mode=1)

	redshift = int(d0['ds']['Redshift']*1000)/1000
	
	#p_index = 2.5
	#lsyn = luminoisty_syn(sn, ncore=ncore,box=[[low]*3,[up]*3],Tsh=Tsh,p_index=p_index)
	#print('Synchrotron normalization factor: {} [z, A, p]'.format(lsyn))
	#lp = np.linspace(2.0,3.0,11)
	#llsyn = np.array([luminosity_syn(sn, ncore=ncore,box=[[low]*3,[up]*3],Tsh=Tsh,p_index=x,facB=1e-4) for x in lp]).T
	#totxt('Lp_syn_wdm_'+str(sn)+'.txt', llsyn,0,0,0)

	nump = 32
	lnu = 10**np.linspace(np.log10(numin())+0.5,np.log10(numax())+1.0,nump)
	np_core = int(nump/ncore)
	lpr = [[i*np_core, (i+1)*np_core] for i in range(ncore-1)] + [[(ncore-1)*np_core, nump]]
	manager = mp.Manager()
	output = manager.Queue()
	def sess(low, up, i):
		output.put([i]+[intensity(d0, lnu[x])['L'] for x in range(low, up)])
	processes = [mp.Process(target=sess, args=(lpr[i][0], lpr[i][1], i)) for i in range(ncore)]
	for p in processes:
		p.start()
	for p in processes:
		p.join()
	raw = [output.get() for p in processes]
	raw = sorted(raw, key=lambda x: x[0])
	lmap = np.hstack([x[1:] for x in raw])
	plt.plot(lnu,lmap)
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel(r'$\nu\ [\mathrm{Hz}]$')
	plt.ylabel(r'$L_{\nu}\ [\mathrm{erg\ s^{-1}\ Hz^{-1}}]$')
	plt.title(lmodel[indm]+r': $z='+str(redshift)+'$',size=14)
	plt.tight_layout()
	plt.savefig('luminosity_'+str(lmodel[indm])+'_'+str(sn)+'_'+str(bins)+'.pdf')
	#plt.show()
	totxt('Lnu_wdm_'+str(sn)+'.txt',[lnu,lmap],0,0,0)
	print('Peak luminosity: {} [erg s^-1 Hz^-1]'.format(max(lmap)))

	nu = 0.1#int(numax()/100/1e9)
	test = intensity(d0, nu*1e9)
	test_ = SZy(d0)
	plt.figure()

	#plt.contourf(test['Y'],test['X'],np.log10(test['I']),np.linspace(-36,-19,100),cmap=plt.cm.gist_ncar)#norm=LogNorm()
	plt.contourf(test['Y'],test['X'],np.log10(test['I']),np.linspace(min(-32.5,np.max(np.log10(test['I']))),max(-16.5,np.max(np.log10(test['I']))),100),cmap=plt.cm.gist_ncar)#norm=LogNorm()	
	cb = plt.colorbar()
	cb.set_label(r'$\log(I_{\nu}\ [\mathrm{erg\ s^{-1}\ cm^{-2}\ Hz^{-1}\ sr^{-1}}])$')
	plt.xlabel(r'$x\ [h^{-1}\mathrm{kpc}]$')
	plt.ylabel(r'$y\ [h^{-1}\mathrm{kpc}]$')
	#plt.title(lmodel[indm]+r': $z='+str(redshift)+'$'+r', $\nu='+str(nu)+'\ \mathrm{GHz}$',size=14)
	plt.tight_layout()
	plt.savefig('intensity_'+str(lmodel[indm])+'_'+str(sn)+'_'+str(bins)+'.pdf')
	#plt.show()

	plt.figure()
	plt.contourf(test_['Y'],test_['X'],np.log10(test_['SZ']),np.linspace(-19.5,max(-9.5,np.max(np.log10(test_['SZ']))),100),cmap=plt.cm.gnuplot)#norm=LogNorm()
	cb = plt.colorbar()
	cb.set_label(r'$\log(y_{\mathrm{LoS}})$')
	plt.xlabel(r'$x\ [h^{-1}\mathrm{kpc}]$')
	plt.ylabel(r'$y\ [h^{-1}\mathrm{kpc}]$')
	plt.title(lmodel[indm]+r': $z='+str(redshift)+'$',size=14)
	plt.tight_layout()
	plt.savefig('SZ_'+str(lmodel[indm])+'_'+str(sn)+'_'+str(bins)+'.pdf')
	#plt.show()

	h = 0.6774
	z = d0['ds']['Redshift']
	dc = DZ(z)*(1+z)
	dz = d0['bins'][2]*UL/(1+z)/h
	angle = (dz/dc)**2
	Vc = (d0['bins'][0] * d0['bins'][1] * d0['bins'][2])*(UL/h/(1+z))**3
	Imax = np.max(test['I'])
	print('Flux: {} [erg s^-1 cm^-2 Hz^-1]'.format((1+z)*angle*np.sum(test['I'])))
	print('Peak intensity: {} [Jy sr^-1]'.format(Imax*1e23/(1+z)**3))
	print('size: {} ["]'.format(np.sum(test['I']>Imax/10)**0.5 *d0['bins'][0]*UL/h/DZ(z) *180/np.pi *60**2))
	#print('Flux1: {} [erg s^-1 cm^-2 Hz^-1]'.format(angle*dz*np.sum(d1['eps'])/4/np.pi ))
	print('Power: {} [erg s^-1]'.format(luminosity(d0)[1]))
	print('Average yLoS: {}'.format(np.average(test_['SZ'])))
	#print('Power1: {} [erg s^-1 Hz^-1]'.format(Vc*np.sum(d1['eps'])))
	
	rne = [-10, 0]
	rT = [int(np.log10(Tsh)), 5]
	reps = [-58, -42]
	histb = 200

	plt.figure()
	a = plt.hist(np.log10(d0['ne'][d0['ne']>1e-10]),bins=histb,alpha=0.5,label='CIC',range=rne)#,density=True)
	#b = plt.hist(np.log10(d1['ne'][d1['ne']>1e-10]),bins=histb,alpha=0.5,label='SPH grid',range=rne)#,density=True)
	plt.xlabel(r'$\log(n_{\mathrm{e}}\ [\mathrm{cm^{-3}}])$')
	#plt.ylabel(r'Probability density')
	dh = (rne[1]-rne[0])/histb
	plt.ylabel(r'$\frac{d N_{\mathrm{grid}}}{d\log n_{\mathrm{e}}}\ [\mathrm{('+str(dh)+'\ dex)^{-1}}]$')
	#plt.yscale('log')
	plt.legend()
	plt.tight_layout()
	plt.savefig('logne_'+str(lmodel[indm])+'_'+str(sn)+'_'+str(bins)+'.pdf')
	#plt.show()
	
	plt.figure()
	a = plt.hist(np.log10(d0['T'][d0['T']>30]),bins=histb,alpha=0.5,label='CIC',range=rT)#,density=True)
	#b = plt.hist(np.log10(d1['T'][d1['T']>30]),bins=histb,alpha=0.5,label='SPH grid',range=rT)#,density=True)
	plt.xlabel(r'$\log(T\ [\mathrm{K}])$')
	#plt.ylabel(r'Probability density')
	dh = (rT[1]-rT[0])/histb
	plt.ylabel(r'$\frac{d N_{\mathrm{grid}}}{d\log T}\ [\mathrm{('+str(dh)+'\ dex)^{-1}}]$')
	#plt.yscale('log')
	plt.legend()
	plt.tight_layout()
	plt.savefig('logT_'+str(lmodel[indm])+'_'+str(sn)+'_'+str(bins)+'.pdf')
	#plt.show()

	plt.figure()
	a = plt.hist(np.log10(d0['eps'][d0['eps']>1e-60]),bins=histb,alpha=0.5,label='CIC',range=reps)#,density=True)
	#b = plt.hist(np.log10(d1['eps'][d1['eps']>1e-60]),bins=histb,alpha=0.5,label='SPH grid',range=reps)#,density=True)
	plt.xlabel(r'$\log(\epsilon_{\nu}\ [\mathrm{erg\ s^{-1}\ cm^{-3}\ Hz^{-1}}])$')
	#plt.ylabel(r'Probability density')
	dh = (reps[1]-reps[0])/histb
	plt.ylabel(r'$\frac{d N_{\mathrm{grid}}}{d\log \epsilon_{\nu}}\ [\mathrm{('+str(dh)+'\ dex)^{-1}}]$')
	#plt.yscale('log')
	#plt.xlim(-50,-40)
	plt.legend()
	plt.tight_layout()
	plt.savefig('logeps_'+str(lmodel[indm])+'_'+str(sn)+'_'+str(bins)+'.pdf')
	#plt.show()
	

