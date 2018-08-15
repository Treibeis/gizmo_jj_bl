from radio import *

if __name__ == "__main__":
	rep0 = 'halo1_jj/'
	lL_nu0 = retxt(rep0+'Lnu_cdm.txt',2,0,0)
	L_nu0 = interp1d(np.log10(lL_nu0[0]),lL_nu0[1])

	lL_nu1 = retxt(rep0+'Lnu_wdm.txt',2,0,0)
	L_nu1 = interp1d(np.log10(lL_nu1[0]),lL_nu1[1])

	lLp0 = retxt(rep0+'Lp_syn_cdm.txt',3,0,0)
	Lp0 = interp1d(lLp0[2],np.log10(lLp0[1]))

	lLp1 = retxt(rep0+'Lp_syn_wdm.txt',3,0,0)
	Lp1 = interp1d(lLp1[2],np.log10(lLp1[1]))

	p_index = 2.198
	A_syn0 = 10**Lp0(p_index)#1.1202483164633895e+59
	A_syn1 = 10**Lp1(p_index)#1.8152757384235942e+58
	facB = 1.0
	facn = 3.4619e-5
	#facn = 1.0
	L_nu_syn0 = lambda x: jnu_syn_(10**x, A_syn0, p_index)*facB**((p_index+1)/4)*facn
	L_nu_syn1 = lambda x: jnu_syn_(10**x, A_syn1, p_index)*facB**((p_index+1)/4)*facn

	redshift = lLp0[0][0]
	dL = DZ(redshift)*(1+redshift)
	print('Syn Flux: {} [erg s^-1 cm^-2 Hz^-1]'.format((1+redshift)*L_nu_syn0(np.log10(1.4e9*(1+redshift)))/dL**2/4/np.pi))

	hmf00 = hmf.MassFunction()
	hmf00.update(n=0.966, sigma_8=0.829,cosmo_params={'Om0':0.315,'H0':67.74},Mmin=3,Mmax=9)
	hmf0 = hmf.MassFunction()
	hmf0.update(n=0.966, sigma_8=0.829,cosmo_params={'Om0':0.315,'H0':67.74},Mmin=8,Mmax=11)
	hmf1 = hmf.wdm.MassFunctionWDM()
	hmf1.update(n=0.966, sigma_8=0.829,cosmo_params={'Om0':0.315,'H0':67.74},Mmin=8,Mmax=11,wdm_mass=3)

	nzbin = 61
	#lt_base=np.linspace(TZ(30),TZ(0),31)/1e9/YR
	lz_base=np.linspace(0,30,nzbin)#np.array([ZT(np.log10(x)) for x in lt_base])
	ln_M_z0 = np.zeros(nzbin)
	ln_M_z1 = np.zeros(nzbin)
	ln_M_z00 = np.zeros(nzbin)
	h = 0.6774
	for i in range(len(lz_base)):
		hmf0.update(z=lz_base[i])
		lm = np.log10(hmf0.m/h)
		ln = np.log10(hmf0.ngtm*h**3)
		nm = interp1d(lm,ln)
		ln_M_z0[i] = (10**nm(9)-10**nm(10))
		hmf1.update(z=lz_base[i])
		lm = np.log10(hmf1.m/h)
		ln = np.log10(hmf1.ngtm*h**3)
		nm = interp1d(lm,ln)
		ln_M_z1[i] = (10**nm(9)-10**nm(10))
		hmf00.update(z=lz_base[i])
		lm = np.log10(hmf00.m/h)
		ln = np.log10(hmf00.ngtm*h**3)
		nm = interp1d(lm,ln)
		ln_M_z00[i] = (10**nm(-2*np.log10((1+lz_base[i])/10)+6)-10**nm(np.log10(2.5*((1+lz_base[i])/10)**-1.5)+7))
	hmf00.update(z=6)
	lm = np.log10(hmf00.m/h)
	ln = np.log10(hmf00.ngtm*h**3)
	nm = interp1d(lm,ln)
	ln_M_z_norm = (10**nm(-2*np.log10((1+6)/10)+6)-10**nm(np.log10(2.5*((1+6)/10)**-1.5)+7))

	lJ21_z = [10.5*((1+lz_base[i])/7)**-1.5*ln_M_z00[i]/ln_M_z_norm for i in range(len(lz_base))]
	J21_z = interp1d(lz_base,np.log10(lJ21_z))
#J21_z = lambda z: J21_z0(z)*(z>6) + np.log10(10**J21_z0(6)*np.exp((1/7-1/(1+z))*20))*(z<=6)
	#"""
	lt_base = [TZ(x)/YR/1e9 for x in lz_base]
	plt.figure()
	plt.plot(lz_base,ln_M_z1,label=lmodel[1])
	plt.plot(lz_base,ln_M_z0,'--',label=lmodel[0])
	plt.xlabel(r'$z$')
	#plt.xlabel(r'$t\ [\mathrm{Gyr}]$')
	plt.ylabel(r'$n(10^{9}-10^{10}\ M_{\odot})\ [h^{3}\mathrm{Mpc^{-3}}]$')
	plt.xlim(0,20)
	plt.ylim(0,3.5)
	plt.legend()
	plt.tight_layout()
	plt.savefig(rep0+'n_M_z.pdf')
	plt.figure()
	plt.plot(lt_base,ln_M_z1,label=lmodel[1])
	plt.plot(lt_base,ln_M_z0,'--',label=lmodel[0])
	#plt.xlabel(r'$z$')
	plt.xlabel(r'$t\ [\mathrm{Gyr}]$')
	plt.ylabel(r'$n(10^{9}-10^{10}\ M_{\odot})\ [h^{3}\mathrm{Mpc^{-3}}]$')
	#plt.xlim(0,20)
	plt.ylim(0,3.5)
	plt.legend()
	plt.tight_layout()
	plt.savefig(rep0+'n_M_t.pdf')
	#plt.show()
	#"""
	def n_a_CDM(a,h=0.6774):
		z = 1/a-1
		nz = interp1d(lz_base,np.log10(ln_M_z0))
		return 10**nz(z)
	
	def n_a_CDM0(a,h=0.6774):
		z = 1/a-1
		nz = interp1d(lz_base,np.log10(ln_M_z00))
		return 10**nz(z)

	def n_a_WDM(a,h=0.6774):
		z = 1/a-1
		nz = interp1d(lz_base,np.log10(ln_M_z1))
		return 10**nz(z)

	JHII_z = lambda z: 0.3*(z<=6) + 0.3*(z>6)*n_a_CDM0(1/(1+z))/n_a_CDM0(1/7)

	mode = int(sys.argv[1])
	if len(sys.argv)>=3:
		tmode = int(sys.argv[2])
	else:
		tmode = 0
	if tmode==0:
		zend = 0.0
	else:
		zend = 7.5

	lnu = 10**np.linspace(np.log10(50),np.log10(1400),100)
	lJnu_syn0 = [Jnu_cosmic(zend,L=L_nu_syn0,nu=x,n_a=n_a_CDM,mode=tmode) for x in lnu]
	lJnu_syn1 = [Jnu_cosmic(zend,L=L_nu_syn1,nu=x,n_a=n_a_CDM,mode=tmode) for x in lnu]
	plt.figure()
	plt.plot(lnu, Tnu(lnu,np.array(lJnu_syn1)),label='Structure formation, '+lmodel[1])
	plt.plot(lnu, Tnu(lnu,np.array(lJnu_syn0)),'--',label='Structure formation, '+lmodel[0])
	plt.plot(lnu, 1e3*Tnu_sky(lnu),':',label='ARCADE 2')
	plt.xlabel(r'$\nu_{\mathrm{obs}}\ [\mathrm{MHz}]$')
	plt.ylabel(r'$\langle\delta T_{\mathrm{syn}}\rangle\ [\mathrm{mK}]$')
	plt.legend()
	plt.xlim(50,1400)
	#plt.title(r'Synchrotron: $\beta_{B}='+str(facB)+r'$, $\beta_{n}='+str(facn*1e5)+r'\times 10^{-5}$, $p='+str(p_index)+'$',size=14)
	plt.xscale('log')
	if mode==0:
		plt.yscale('log')
	plt.tight_layout()
	if mode==0:
		plt.savefig(rep0+'logTnu_syn.pdf')
	else:
		plt.savefig(rep0+'Tnu_syn.pdf')
	print('Tnu_syn at 310 MHz: {} mK (CDM)'.format(Tnu(310,Jnu_cosmic(zend,L=L_nu_syn0,nu=310,n_a=n_a_CDM,mode=tmode))))

	lp = np.linspace(2.0,3.0,21)
	lT_syn0 = np.array([Tnu(310,Jnu_cosmic(zend,L=lambda x: jnu_syn_(10**x, 10**Lp0(y), y),nu=310,n_a=n_a_CDM,mode=tmode)) for y in lp])
	lT_syn1 = np.array([Tnu(310,Jnu_cosmic(zend,L=lambda x: jnu_syn_(10**x, 10**Lp1(y), y),nu=310,n_a=n_a_WDM,mode=tmode)) for y in lp])
	plt.figure()
	plt.plot(lp,1e3*Tnu_sky(310)/lT_syn1,label=lmodel[1])
	plt.plot(lp,1e3*Tnu_sky(310)/lT_syn0,'--',label=lmodel[0])
	plt.legend()
	plt.xlim(2,3)
	plt.xlabel(r'$p=2\alpha-3$')
	plt.ylabel(r'$\beta_{B,-4}^{(p+1)/4}\beta_{n} \cdot T_{\mathrm{sky}}/\langle\delta T_{\mathrm{syn}}\rangle(310\ \mathrm{MHz})$')
	if mode==0:
		plt.yscale('log')
	plt.tight_layout()
	if mode==0:
		plt.savefig(rep0+'logRat_syn.pdf')
	else:
		plt.savefig(rep0+'Rat_syn.pdf')
	p0 = 2.198
	print('Ratio (CDM): {}'.format(1e3*Tnu_sky(310)/Tnu(310,Jnu_cosmic(zend,L=lambda x: jnu_syn_(10**x, 10**Lp0(p0), p0),nu=310,n_a=n_a_CDM,mode=tmode))))
	print('Ratio (WDM): {}'.format(1e3*Tnu_sky(310)/Tnu(310,Jnu_cosmic(zend,L=lambda x: jnu_syn_(10**x, 10**Lp1(p0), p0),nu=310,n_a=n_a_WDM,mode=tmode))))

	plt.figure()
	plt.plot(lL_nu1[0], lL_nu1[1], label=lmodel[1])
	plt.plot(lL_nu0[0], lL_nu0[1], '--', label=lmodel[0])
	plt.legend()
	plt.xlabel(r'$\nu\ [\mathrm{Hz}]$')
	plt.ylabel(r'$L_{\nu}\ [\mathrm{erg\ s^{-1}\ Hz^{-1}}]$')
	plt.xscale('log')
	if mode==0:
		plt.yscale('log')
	plt.tight_layout()
	if mode==0:
		plt.savefig(rep0+'logLnu_com.pdf')
	else:
		plt.savefig(rep0+'Lnu_com.pdf')

	lz = np.linspace(19.9,1/(1-5e-2)-1,100)
	lJz0 = [Jnu_cosmic(zend,L=L_nu0,nu=1420/(1+z),n_a=n_a_CDM,mode=tmode) for z in lz]
	lJz1 = [Jnu_cosmic(zend,L=L_nu1,nu=1420/(1+z),n_a=n_a_WDM,mode=tmode) for z in lz]
	lJ0 = [Jnu_cosmic(z,L=L_nu0,nu=310,n_a=n_a_CDM,mode=tmode) for z in lz]
	lJ1 = [Jnu_cosmic(z,L=L_nu1,nu=310,n_a=n_a_WDM,mode=tmode) for z in lz]
	if tmode==1:
		nmax1 = np.max(ln_M_z1)
		zmax1 = lz_base[[i for i in range(len(lz_base)) if ln_M_z1[i]==nmax1][0]]
		nmax0 = np.max(ln_M_z0)
		zmax0 = lz_base[[i for i in range(len(lz_base)) if ln_M_z0[i]==nmax0][0]]
		print('zmax: {} (CDM), {} (WDM)'.format(zmax0, zmax1))
		J1max = Jnu_cosmic(zmax1,L=L_nu1,nu=310,n_a=n_a_WDM,mode=tmode)
		J0max = Jnu_cosmic(zmax0,L=L_nu0,nu=310,n_a=n_a_CDM,mode=tmode)
		for i in range(len(lJ0)):
			if lz[i]<=zmax1:
				lJ1[i] = J1max
			if lz[i]<=zmax0:
				lJ0[i] = J0max

	lnu = 10**np.linspace(np.log10(50),np.log10(1400),100)
	lJnu0 = [Jnu_cosmic(zend,L=L_nu0,nu=x,n_a=n_a_CDM,mode=tmode) for x in lnu]
	lJnu1 = [Jnu_cosmic(zend,L=L_nu1,nu=x,n_a=n_a_WDM,mode=tmode) for x in lnu]
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	#ax2 = ax1.twiny()
	ax1.plot(lnu, Tnu(lnu,np.array(lJnu1)), label=r'Structure formation, '+lmodel[1],lw=1)
	ax1.plot(lnu, Tnu(lnu,np.array(lJnu0)), label=r'Structure formation, '+lmodel[0],ls='--',lw=1)
	ax1.plot(lnu, Tnu(lnu,JHII_z(6)),'-.',label=r'Mini-halo $\mathrm{H_{II}}$ regions',lw=1)# ($M_{*}\sim 100\ M_{\odot})$',lw=1)
	#ax1.plot(lnu[lnu<1420/7], Tnu(lnu[lnu<1420/7],10**J21_z(1420/lnu[lnu<1420/7]-1)),color='r',ls=':',lw=1)
	ax1.plot(lnu[lnu>0], Tnu(lnu[lnu>0],10**J21_z(1420/lnu[lnu>0]-1)), ls=':',color='r', label=r'21 cm emission',lw=1)#, $\nu_{\mathrm{obs}}=1420/(1+z)\ \mathrm{MHz}$',lw=1)
	lTnu_IGM = [Tnu(x,Jnu_bg(x)) for x in lnu]
	ax1.plot(lnu, lTnu_IGM, ls='-.',lw=2,color='g',label=r'ionized diffuse IGM')
	ax1.plot(lnu,Tnu_SKA(lnu),'k--',label=r'SKA',lw=2)#, 10$\sigma$, $10^{3}$ h',lw=2)
	ax1.fill_between(lnu,1e3*Tnu_sky_ff(lnu,-1),1e3*Tnu_sky_ff(lnu,1),facecolor='gray',label=r'$T_{\mathrm{ff}}^{\mathrm{G}}$',alpha=0.5)
	#ax2.set_xscale('log')
	#loc = [1420/3.0,1420/7.0,1420/10.215,1420/13.593,1420/21]
	#ax2.set_xticks(loc)
	#ax2.set_xticklabels(['Post-reionization','6.0','9.2','13.6','20'],size=11)
	#ax2.set_xlabel(r'$z$')#=1420/(1+z)\ \mathrm{MHz}
	yup = np.max([160.0,np.max(Tnu(lnu,np.array(lJnu1))),np.max(Tnu(lnu,np.array(lJnu0)))])*1.05
	ax1.plot([1420/7,1420/7],[1e-3,yup],lw=0.5,color='k')
	#ax1.plot([1420/10.215,1420/10.215],[1e-3,210.0],'--',lw=0.5,color='k')
	#ax1.plot([1420/13.593,1420/13.593],[1e-3,210.0],'--',lw=0.5,color='k')
	ax1.fill_between([1420/7,1400],[1e-3,1e-3],[yup,yup],facecolor='gray',alpha=0.2)
	ax1.set_xlim(50,1400)
	#ax2.set_xlim(ax1.get_xlim())
	ax1.set_ylim(1e-4,yup)
	ax1.set_xlabel(r'$\nu_{\mathrm{obs}}\ [\mathrm{MHz}]$')
	ax1.set_ylabel(r'$\langle\delta T\rangle\ [\mathrm{mK}]$')
	ax1.legend()
	ax1.set_xscale('log')
	if mode==0:
		ax1.set_yscale('log')
	plt.tight_layout()
	if mode==0:
		plt.savefig(rep0+'logTnu.pdf')
	else:
		plt.savefig(rep0+'Tnu.pdf')
	print('Tnu at 310 MHz: {} (CDM) mK'.format(Tnu(310,Jnu_cosmic(zend,L=L_nu0,nu=310,n_a=n_a_CDM,mode=tmode))))
	print('Tnu at 310 MHz: {} (WDM) mK'.format(Tnu(310,Jnu_cosmic(zend,L=L_nu1,nu=310,n_a=n_a_WDM,mode=tmode))))
	#plt.show()

	lnu = 10**np.linspace(np.log10(1e3),np.log10(9e4),100)
	lJnu0 = [Jnu_cosmic(zend,L=L_nu0,nu=x,n_a=n_a_CDM,mode=tmode) for x in lnu]
	lJnu1 = [Jnu_cosmic(zend,L=L_nu1,nu=x,n_a=n_a_WDM,mode=tmode) for x in lnu]
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax1.plot(lnu/1e3, Tnu(lnu,np.array(lJnu1)), label=r'Structure formation, '+lmodel[1],lw=1)
	ax1.plot(lnu/1e3, Tnu(lnu,np.array(lJnu0)), label=r'Structure formation, '+lmodel[0],ls='--',lw=1)
	ax1.plot(lnu/1e3, Tnu(lnu,JHII_z(6)),'-.',label=r'Mini-halo $\mathrm{H_{II}}$ regions',lw=1)# ($M_{*}\sim 100\ M_{\odot})$',lw=1)
	ax1.plot(lnu/1e3, [Tnu(x,Jnu_bg(x)) for x in lnu], ls='-.',lw=2,color='g',label=r'ionized diffuse IGM')
	ax1.plot(lnu/1e3,Tnu_SKA(lnu),'k--',label=r'SKA', lw=2)#, 10$\sigma$, $10^{3}$ h',lw=2)
	ax1.fill_between(lnu/1e3,1e3*Tnu_sky_ff(lnu,-1),1e3*Tnu_sky_ff(lnu,1),facecolor='gray',label=r'$T_{\mathrm{ff}}^{\mathrm{G}}$',alpha=0.5)
	ax1.set_xlim(1.0,90)
	yup = np.max([0.31,np.max(Tnu(lnu,np.array(lJnu1))),np.max(Tnu(lnu,np.array(lJnu0)))])*1.05
	ax1.set_ylim(1e-6,yup)
	ax1.set_xlabel(r'$\nu_{\mathrm{obs}}\ [\mathrm{GHz}]$')
	ax1.set_ylabel(r'$\langle\delta T\rangle\ [\mathrm{mK}]$')
	ax1.legend()
	ax1.set_xscale('log')
	if mode==0:
		ax1.set_yscale('log')
	plt.tight_layout()
	if mode==0:
		plt.savefig(rep0+'logTnu0.pdf')
	else:
		plt.savefig(rep0+'Tnu0.pdf')
	
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	#ax2 = ax1.twiny()
	ax1.plot(lz, Tnu(310,np.array(lJ1)), label=r'$>z$, $\nu_{\mathrm{obs}}=310\ \mathrm{MHz}$, '+lmodel[1],lw=1)
	ax1.plot(lz, Tnu(310,np.array(lJ0)), '--',label=r'$>z$, $\nu_{\mathrm{obs}}=310\ \mathrm{MHz}$, '+lmodel[0],lw=1)
	#ax1.plot(lz, lJ2, '-.', label=r'$>z$, $\nu_{\mathrm{obs}}=10^{4}\ \mathrm{MHz}$')
	#ax1.plot(lz, Tnu(1420/(1+lz),JHII_z(6)),'-.',label=r'$\mathrm{H_{II}}$ regions ($M_{*}\sim 100\ M_{\odot})$',lw=2)
	#ax1.plot(lz[lz>6], Tnu(1420/(1+lz[lz>6]),10**J21_z(lz[lz>6])),lw=2,color='r',ls=':')
	#ax1.plot(lz, Tnu(1420/(1+lz),np.array(lJz1)), label=r'overall ($z='+str(zend)+'$), '+lmodel[1],lw=2,color='k')
	#ax1.plot(lz, Tnu(1420/(1+lz),np.array(lJz0)), label=r'overall ($z='+str(zend)+'$), '+lmodel[0],lw=2,color='k',ls='--')
	#ax1.plot(lz[lz>=0], Tnu(1420/(1+lz[lz>=0]),10**J21_z(lz[lz>=0])), ls=':',lw=2,color='r', label=r'21-cm emission, $\nu_{\mathrm{obs}}=1420/(1+z)\ \mathrm{MHz}$')
	#ax1.plot(lz, [Tnu(1420/(1+x),Jnu_bg(1420/(1+x))) for x in lz], ls='-.',lw=2,color='g',label=r'ionized diffuse IGM')
	loc0 = np.linspace(7.5,20.0,6)
	loc = np.hstack([[2.0,6.0],loc0])
	#ax2.set_xticks(loc)
	#ax2.set_xticklabels(['Post-reionization','203']+[str(int(x)) for x in 1420/(loc0+1)],size=11)
	#ax2.set_xlabel(r'$\nu_{\mathrm{obs}}\ [\mathrm{MHz}]$')#=1420/(1+z)\ \mathrm{MHz}
	#print([min(lJ2),max(lJ1)])
	yup = np.min([6e3, np.max([60, np.max(Tnu(1420/(1+lz),np.array(lJz1))),np.max(Tnu(1420/(1+lz),np.array(lJz0))), np.max(Tnu(310,np.array(lJ0)))])*1.05])
	ax1.plot([6,6],[1e-9,yup],lw=0.5,color='k')
	#ax1.plot([9.215,9.215],[1e-9,60],'--',lw=0.5,color='k')
	#ax1.plot([12.593,12.593],[1e-9,60],'--',lw=0.5,color='k')
	ax1.fill_between([0,6],[1e-9,1e-9],[yup,yup],facecolor='gray',alpha=0.2)
	#ax2.set_title(r'$\langle L_{\nu}\rangle \sim 1-4\times 10^{24}\ \mathrm{erg\ s^{-1}\ Hz^{-1}}$, $10^{9}\lesssim M_{\mathrm{halo}}\ [h^{-1}M_{\odot}]\lesssim 10^{10}$,'+'\n'+r' $\nu_{\mathrm{obs}}\sim 1- 10^{7}\ [\mathrm{MHz}/(1+z)]$',size=12)#$10/(1+z)\lesssim \nu_{\mathrm{obs}}\ [\mathrm{MHz}]\lesssim 10^{7}/(1+z)$')
	ax1.set_xlim(0,20)
	#ax2.set_xlim(ax1.get_xlim())
	ax1.set_ylim(1e-9,yup)
	ax1.set_xlabel(r'$z_{\mathrm{end}}$')
	ax1.set_ylabel(r'$\langle\delta T\rangle(>z_{\mathrm{end}}) [\mathrm{mK}]$')
	ax1.legend()
	if mode==0:
		ax1.set_yscale('log')
	plt.tight_layout()
	if mode==0:
		plt.savefig(rep0+'logTnuz.pdf')
	else:
		plt.savefig(rep0+'Tnuz.pdf')
	lT_ = Tnu(310,np.array(lJ0))
	zmax = lz[[i for i in range(len(lT_)) if lT_[i]==np.max(lT_)][0]]
	print('Tnu at 310 MHz: {}, for z = {}'.format(np.max(lT_), zmax))
	#plt.show()

	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax2 = ax1.twiny()
	ax1.plot(lz, lJ1, label=r'$>z$, $\nu_{\mathrm{obs}}=310\ \mathrm{MHz}$, '+lmodel[1],lw=1)
	ax1.plot(lz, lJ0, '--',label=r'$>z$, $\nu_{\mathrm{obs}}=310\ \mathrm{MHz}$, '+lmodel[0],lw=1)
	#ax1.plot(lz, lJ2, '-.', label=r'$>z$, $\nu_{\mathrm{obs}}=10^{4}\ \mathrm{MHz}$')
	ax1.plot(lz, JHII_z(lz),'-.',label=r'$>z$, $\mathrm{H_{II}}$ regions ($M_{*}\sim 100\ M_{\odot})$',lw=1)
	#ax1.plot(lz[lz>6], 10**J21_z(lz[lz>6]),lw=2,color='r',ls=':')
	ax1.plot(lz, lJz1, label=r'overall ($z='+str(zend)+'$), '+lmodel[1],lw=2,color='k')
	ax1.plot(lz, lJz0, label=r'overall ($z='+str(zend)+'$), '+lmodel[0],lw=2,color='k',ls='--')
	ax1.plot(lz[lz>=0], 10**J21_z(lz[lz>=0]), ls=':',lw=2,color='r', label=r'21-cm emission, $\nu_{\mathrm{obs}}=1420/(1+z)\ \mathrm{MHz}$')
	ax1.plot(lz, [Jnu_bg(1420/(1+x)) for x in lz], ls='-.',lw=2,color='g',label=r'ionized diffuse IGM')
	loc0 = np.linspace(7.5,20.0,6)
	loc = np.hstack([[2.0,6.0],loc0])
	ax2.set_xticks(loc)
	ax2.set_xticklabels(['Post-reionization','203']+[str(int(x)) for x in 1420/(loc0+1)],size=11)
	ax2.set_xlabel(r'$\nu_{\mathrm{obs}}\ [\mathrm{MHz}]$')#=1420/(1+z)\ \mathrm{MHz}
	#print([min(lJ2),max(lJ1)])
	yup = np.max([22.5, np.max(lJz1), np.max(lJz0), np.max(lJ0)])*1.05
	ax1.plot([6,6],[1e-9,yup],lw=0.5,color='k')
	#ax1.plot([9.215,9.215],[1e-9,22.5],'--',lw=0.5,color='k')
	#ax1.plot([12.593,12.593],[1e-9,22.5],'--',lw=0.5,color='k')
	ax1.fill_between([0,6],[1e-9,1e-9],[yup,yup],facecolor='gray',alpha=0.2)
	#plt.plot([0.42,0.42],[1e-9,20],ls='--',lw=0.5,color='k')
	#ax2.set_title(r'$\langle L_{\nu}\rangle \sim 1-4\times 10^{24}\ \mathrm{erg\ s^{-1}\ Hz^{-1}}$, $10^{9}\lesssim M_{\mathrm{halo}}\ [h^{-1}M_{\odot}]\lesssim 10^{10}$,'+'\n'+r' $\nu_{\mathrm{obs}}\sim 1- 10^{7}\ [\mathrm{MHz}/(1+z)]$',size=12)#$10/(1+z)\lesssim \nu_{\mathrm{obs}}\ [\mathrm{MHz}]\lesssim 10^{7}/(1+z)$')
	ax1.set_xlim(0,20)
	ax2.set_xlim(ax1.get_xlim())
	ax1.set_ylim(1e-9,yup)
	ax1.set_xlabel(r'$z$')
	ax1.set_ylabel(r'$J_{\nu_{\mathrm{obs}}}(>z)\ [\mathrm{Jy\ sr^{-1}}]$')
	ax1.legend()
	if mode==0:
		ax1.set_yscale('log')
	plt.tight_layout()
	if mode==0:
		plt.savefig(rep0+'logJnu.pdf')
	else:
		plt.savefig(rep0+'Jnu.pdf')
	
	lz = np.linspace(20.0,0,100)
	rat0 = np.array([Jnu_cosmic(zend,L=L_nu0,nu=1420/(1+z),n_a=n_a_CDM,mode=tmode)/10**J21_z(z) for z in lz])
	rat1 = np.array([Jnu_cosmic(zend,L=L_nu1,nu=1420/(1+z),n_a=n_a_WDM,mode=tmode)/10**J21_z(z) for z in lz])
	lJ3 = [Jnu_cosmic(z,L=L_nu1,nu=80,n_a=n_a_WDM,mode=tmode) for z in lz]
	#lJ4 = [Jnu_cosmic(z,L=L_nu1,nu=139,n_a=n_a_WDM) for z in lz]
	lJ5 = [Jnu_cosmic(z,L=L_nu0,nu=80,n_a=n_a_CDM,mode=tmode) for z in lz]
	#lJ6 = [Jnu_cosmic(z,L=L_nu0,nu=104.5,n_a=n_a_CDM) for z in lz]
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax2 = ax1.twiny()
	loc0 = np.linspace(7.5,20.0,6)
	loc = np.hstack([[2.0,6.0],loc0])
	ax2.set_xticks(loc)
	ax2.set_xticklabels(['Post-reionization','203']+[str(int(x)) for x in 1420/(loc0+1)],size=11)
	ax2.set_xlabel(r'$\nu_{\mathrm{obs}}\ [\mathrm{MHz}]$')#=1420/(1+z)\ \mathrm{MHz}
	ax1.plot(lz, lJ3/10**J21_z(1420/80-1),label=r'$>z$, $\nu_{\mathrm{obs}}=80\ \mathrm{MHz}$, '+lmodel[1],lw=1)
	ax1.plot(lz, lJ5/10**J21_z(1420/80-1),'--',label=r'$>z$, $\nu_{\mathrm{obs}}=80\ \mathrm{MHz}$, '+lmodel[0],lw=1)
	#ax1.plot(lz, lJ4/10**J21_z(1420/139-1),'-.',label=r'$>z$, $\nu_{\mathrm{obs}}=139\ \mathrm{MHz}$, '+lmodel[1],lw=1)
	#ax1.plot(lz, lJ6/10**J21_z(1420/104.5-1),':',label=r'$>z$, $\nu_{\mathrm{obs}}=104.5\ \mathrm{MHz}$, '+lmodel[0],lw=1)
	ax1.plot(lz[lz>6], rat1[lz>6],'k-',label=r'overall ($z='+str(zend)+'$), '+lmodel[1],lw=2)#, $\nu_{\mathrm{obs}}=1420/(1+z)\ \mathrm{MHz}$
	ax1.plot(lz[lz<6], rat1[lz<6],'k:',lw=2)
	ax1.plot(lz[lz>6], rat0[lz>6],'k--',label=r'overall ($z='+str(zend)+'$), '+lmodel[0],lw=2)#, $\nu_{\mathrm{obs}}=1420/(1+z)\ \mathrm{MHz}$
	ax1.plot(lz[lz<6], rat0[lz<6],'k:',lw=2)
	ax1.plot(lz, [1 for x in lz],'k:',lw=0.5)
	yup = np.max([2.5, np.max(rat1[lz>6]), np.max(rat0[lz>6]), np.max(lJ5/10**J21_z(1420/80-1))])*1.05
	ax1.plot([6,6],[0,yup],lw=0.5,color='k')
	#ax1.plot([9.215,9.215],[0,12],'--',lw=0.5,color='k')
	#ax1.plot([12.593,12.593],[0,12],'--',lw=0.5,color='k')
	ax1.fill_between([0,6],[0,0],[yup,yup],facecolor='gray',alpha=0.2)
	ax1.set_xlim(0,20)
	ax1.set_ylim(1e-2,yup)
	#plt.yscale('log')
	ax1.set_xlabel(r'$z$')
	ax1.set_ylabel(r'$\frac{J_{\nu_{\mathrm{obs}},\mathrm{Structure\ formation}}(>z)}{J_{\nu_{\mathrm{obs}},\mathrm{21\ cm}}}$')
	ax1.legend()
	ax1.set_yscale('log')
	plt.tight_layout()
	plt.savefig(rep0+'Jnu_ratio.pdf')
	#plt.show()

	

