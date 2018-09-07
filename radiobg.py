from radio import *
d_delta = lambda z: 1.686*(1-0.01*(1+z)/20)
h = 0.6774

Mmax = 10
Mref = 7e9
NUREF = 1e11
BETA_l = 0.0#5/3
BETA_t = 0.0#1-5/3
Tref = (TZ(7.22147651692)-TZ(10.2385321471))/YR

hmf000 = hmf.MassFunction()
hmf000.update(n=0.966, sigma_8=0.829,cosmo_params={'Om0':0.315,'H0':67.74},Mmin=3,Mmax=9)
hmf000.update(z=6)
lm = np.log10(hmf000.m/h)
ln = np.log10(hmf000.ngtm)#*h**3)
nm = interp1d(lm,ln)
Nps_ref = (10**nm(-2*np.log10((1+6)/10)+6)-10**nm(np.log10(2.5*((1+6)/10)**-1.5)+7))/(UL*1e3/h)**3
Nion = 3e-21*2.6e-13/(SPEEDOFLIGHT*1e-39*Nps_ref)

def extinction(z, m, nu, delta = 200):
	T = Tvir(m, z)
	M = m*UM/1e10
	R = (M/(rhom(1/(1+z))*delta)*3/4/np.pi)**(1/3)
	od = 4*CHARGE**6/(3*ELECTRON*SPEEDOFLIGHT*BOL) *(2*np.pi/(3*BOL*ELECTRON))**0.5 *T**-1.5 * (delta*rhom(1/(1+z))/(PROTON*mmw()))**2 * (nu*1e6)**-2 * R
	return od

def tau_mini(z):
	return 1/(2.6e-13*rhom(1/(1+z))/(PROTON*mmw()))/YR

Mup = lambda z: 2.5e7*((1+z)/10)**-1.5
Mdown = lambda z: 1e6*((1+z)/10)**-2
Mmini = lambda z: (Mup(z)+Mdown(z))/2

def Nion_m(m):
	return Nion*m/Mmini(6)

def tau_M(z, Mref = Mref, tref = Tref, Mbd = 1e10, beta = BETA_t):
	alpha = np.log(tref/tau_mini(z))/np.log(Mref/Mup(z))
	def func(m):
		if m<=Mup(z):
			return tau_mini(z)
		elif m<=Mbd:
			return tref*(m/Mref)**alpha
		else:
			return tref*(Mbd/Mref)**alpha *(m/Mbd)**beta
		#return 100e6
	return func

def Lnu_minih(m, z, Tmini = 1e3):
	return Nion_m(m)*rhom(1/(1+z))/(PROTON*mmw()) * Tmini**-0.5 * 2**5*np.pi*CHARGE**6/3/ELECTRON/SPEEDOFLIGHT**3 * (2*np.pi/3/BOL/ELECTRON)**0.5 

def Lnu_M(L, z, Mref = Mref, Mbd = 1e10, lognu_ref = np.log10(NUREF), Tmini = 1e3, beta = BETA_l, delta=200, nu_min = 1e3):
	mup = Mup(z)
	alpha = np.log(L(lognu_ref)/Lnu_minih(mup,z))/np.log(Mref/mup)
	#print(alpha)
	def func(m, nu):
		if m<=mup:
			return Lnu_minih(m,z) * np.exp(-HBAR*2*np.pi*nu*1e6/BOL/Tmini)
		elif m<=Mbd:
			return L(np.log10(nu*1e6))*(m/Mref)**alpha * np.exp(-HBAR*2*np.pi*nu*1e6/BOL/Tvir(m, z)) #* np.exp(-extinction(z, m, nu))
		else:
			R = (m/(rhom(1/(1+z))*delta)*3/4/np.pi)**(1/3)
			Rref = (Mref/(rhom(1/(1+z))*delta)*3/4/np.pi)**(1/3)
			nu_bd = nu_min * (R/Rref)**0.5
			if nu<nu_bd:
				return L(np.log10(nu_bd*1e6))*(Mbd/Mref)**alpha * (m/Mbd)**beta * (nu/nu_bd)**2
			else:
				return L(np.log10(nu*1e6))*(Mbd/Mref)**alpha * (m/Mbd)**beta #* np.exp(-HBAR*2*np.pi*nu*1e6/BOL/Tvir(m, z))
	return func

from scipy.interpolate import interp2d

def dndm_z(z1=0, z0=31, mode=0, nbin=100, Mmax=10, load=0):
	if load==0:
		lz = np.linspace(z1, z0, nbin)
		out = []
		if mode==0:
			hmf_ = hmf.MassFunction()
			hmf_.update(n=0.966, sigma_8=0.829,cosmo_params={'Om0':0.315,'H0':67.74},Mmin=np.log10(Mdown(max(z0,z1)))-1,Mmax=Mmax+1)
		else:
			hmf_ = hmf.wdm.MassFunctionWDM(wdm_mass=3)
			hmf_.update(n=0.966, sigma_8=0.829,cosmo_params={'Om0':0.315,'H0':67.74},Mmin=np.log10(Mdown(max(z0,z1)))-1,Mmax=Mmax+1)
		for z in lz:
			hmf_.update(z=z)
			out.append(hmf_.dndlog10m)
		lm = hmf_.m
		totxt('mlist.txt',[lm],0,0,0)
		totxt('zlist.txt',[lz],0,0,0)
		totxt('dndm_'+lmodel[mode]+'.txt',out,0,0,0)
	else:
		lm = retxt('mlist.txt',1,0,0)[0]
		lz = retxt('zlist.txt',1,0,0)[0]
		out = retxt('dndm_'+lmodel[mode]+'.txt',nbin,0,0)
	return interp2d(lm,lz,out)

dndm0 = dndm_z(mode=0,Mmax=Mmax,load=1)
dndm1 = dndm_z(mode=1,Mmax=Mmax,load=1)

def Jnu_final(z1, nu, z0 = 30, L = lambda x:1e30, Mref=Mref, dndm=dndm0, zstep = 1.0, h=0.6774):
	start = time.time()
	nzb = int(abs(z0-z1)/zstep)+1
	unit = YR*SPEEDOFLIGHT /(UL*1e3/h)**3/4/np.pi
	lz = np.linspace(z1,z0,nzb)
	ljnu = []
	for z in lz:
		mup = Mup(z)
		mdown = Mdown(z)
		lnu_m = Lnu_M(L, z, Mref)#Lnu_M(L, Lnu_minih(z), Mref, (mup+mdown)/2, z)
		tau = tau_M(z, Mref)#tau_M(200e6, tau_mini(z), Mref, (mup+mdown)/2)
		def dndm_(z, m):
			return dndm(m, z)
		def integrand(m):
			return lnu_m(10**m,nu*(1+z)) * tau(10**m) * max(0.0,-derivative(dndm_, z, 1e-2, args=(10**m,)))
		ljnu.append(quad(integrand, np.log10(mdown), Mmax, epsrel=-4)[0])
	jnu_z = interp1d(lz, ljnu)
	out = quad(jnu_z, z1, z0, epsrel=-4)[0]*unit*1e23
	print('zend = {}, nu = {} [MHz], Jnu = {} [Jy], time: {} s'.format(z1, nu, out, time.time()-start))
	return out

def meanL(a=5/3, b=2.5, g=0.0, m=7):
	return 1#m**(-a+g) * (1-b)*(10**(a-b-g+1)-1) /((a-b-g+1)*(10**(1-b)-1))

if __name__ == "__main__":
	load = 1
	tag = 1
	nbin = 50
	sn_min = 16
	sn_max = 25
	rep0 = 'halo1_jj/'

	mode = int(sys.argv[1])
	if len(sys.argv)>=3:
		zend = float(sys.argv[2])
	else:
		zend = 7.5

	lnu_z0 = np.array(retxt(rep0+'luminosity_z_100_CDM.txt',2,1,0))
	lt0 = np.array([TZ(x)/YR/1e6 for x in lnu_z0[0][sn_min-1:sn_max+1]])
	ldt0 = np.abs(lt0[1:]-lt0[:-1])
	lL_nu0 = retxt(rep0+'Lnu_cdm.txt',2,0,0)
	nbin_nu = len(lL_nu0[0])
	lnu_raw0 = np.zeros(nbin_nu)
	for i in range(sn_min, sn_max+1):
		lnu_raw0 += np.array(retxt(rep0+'Lnu_cdm_'+str(i)+'.txt',2,0,0)[1])*ldt0[i-sn_min]
	lnu0 = lnu_raw0/np.sum(ldt0)

	lnu_z1 = np.array(retxt(rep0+'luminosity_z_100_WDM_3_kev.txt',2,1,0))
	lt1 = np.array([TZ(x)/YR/1e6 for x in lnu_z1[0][sn_min-1:]])
	ldt1 = np.abs(lt1[1:]-lt1[:-1])
	lL_nu1 = retxt(rep0+'Lnu_wdm.txt',2,0,0)
	nbin_nu = len(lL_nu1[0])
	lnu_raw1 = np.zeros(nbin_nu)
	for i in range(sn_min, sn_max+1):
		lnu_raw1 += np.array(retxt(rep0+'Lnu_wdm_'+str(i)+'.txt',2,0,0)[1])*ldt1[i-sn_min]
	lnu1 = lnu_raw1/np.sum(ldt1)

	L_nu0 = interp1d(np.log10(lL_nu0[0]),np.array(lnu0))
	L_nu1 = interp1d(np.log10(lL_nu1[0]),np.array(lnu1))

	print(Tref, np.sum(ldt1))

	plt.figure()
	plt.plot(lL_nu1[0], lL_nu1[1], label=lmodel_[1])
	plt.plot(lL_nu0[0], lL_nu0[1], '--', label=lmodel_[0])
	plt.legend()
	plt.xlabel(r'$\nu\ [\mathrm{Hz}]$')
	plt.ylabel(r'$L^{\mathrm{S}}_{\nu}\ [\mathrm{erg\ s^{-1}\ Hz^{-1}}]$')
	plt.xscale('log')
	if mode==0:
		plt.yscale('log')
	plt.tight_layout()
	if mode==0:
		plt.savefig(rep0+'logLnu_com.pdf')
	else:
		plt.savefig(rep0+'Lnu_com.pdf')

	plt.figure()
	plt.plot(lL_nu1[0], lnu0, label=lmodel_[1])
	plt.plot(lL_nu0[0], lnu1, '--', label=lmodel_[0])
	plt.legend()
	plt.xlabel(r'$\nu\ [\mathrm{Hz}]$')
	plt.ylabel(r'$L^{\mathrm{ref}}_{\nu}\ [\mathrm{erg\ s^{-1}\ Hz^{-1}}]$')
	plt.xscale('log')
	if mode==0:
		plt.yscale('log')
	plt.tight_layout()
	if mode==0:
		plt.savefig(rep0+'logLnu_ref.pdf')
	else:
		plt.savefig(rep0+'Lnu_ref.pdf')

	z_eg = 7.5
	lnu_m0 = Lnu_M(L_nu0, z_eg, 1e10)
	lnu_m1 = Lnu_M(L_nu1, z_eg, 1e10)
	lm = 10**np.linspace(np.log10(Mdown(z_eg)),Mmax, 100)
	lL0 = [lnu_m0(x, NUREF/1e6) for x in lm]
	lL1 = [lnu_m1(x, NUREF/1e6) for x in lm]
	plt.figure()
	plt.plot(lm, lL1, label=lmodel_[0])
	plt.plot(lm, lL0, label=lmodel_[1],ls='--')
	plt.plot(lm, Lnu_minih(lm, z=z_eg),label=r'$L_{\nu}^{\mathrm{mini}}$',ls=':')
	#plt.plot([Mdown(z_eg),Mup(z_eg)], [Lnu_minih(Mdown(z_eg), z_eg),Lnu_minih(Mup(z_eg), z_eg)],marker='*',label='Minihalo',ls=':')
	plt.scatter([Mref],[lnu_m1(Mref, NUREF/1e6)],marker='^',label=r'$L_{\nu}^{\mathrm{ref}}$, '+lmodel_[0],alpha=0.5)
	plt.scatter([Mref],[lnu_m0(Mref, NUREF/1e6)],marker='o',label=r'$L_{\nu}^{\mathrm{ref}}$, '+lmodel_[1],alpha=0.5)
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel(r'$M\ [\odot]$')
	plt.ylabel(r'$L_{\nu='+str(NUREF/1e9)+'\ \mathrm{GHz}}\ [\mathrm{erg\ s^{-1}\ Hz^{-1}}]$')
	plt.legend()
	plt.tight_layout()
	plt.savefig(rep0+'Lnu_M.pdf')
	print('Lnu_ref(CDM): {}, Lnu_ref(WDM): {} [erg s^-1 Hz^-1]'.format(lnu_m0(Mref, NUREF/1e6), lnu_m1(Mref, NUREF/1e6)))
	print('Lnu_mini: {} [erg s^-1 Hz^-1]'.format(Lnu_minih(Mref, z=z_eg)))
	#plt.show()

	lt_m = tau_M(z_eg)
	ltrec = [lt_m(x)/1e6 for x in lm]
	plt.figure()
	plt.plot(lm, ltrec)
	plt.xscale('log')
	#plt.yscale('log')
	plt.xlabel(r'$M\ [\odot]$')
	plt.ylabel(r'$t_{*}\ [\mathrm{Myr}]$')
	plt.tight_layout()
	plt.savefig(rep0+'trec_M.pdf')
	#plt.show()

	if load==0:
		hmf00 = hmf.MassFunction()
		hmf00.update(n=0.966, sigma_8=0.829,cosmo_params={'Om0':0.315,'H0':67.74},Mmin=3,Mmax=9)
		#hmf0 = hmf.MassFunction()
		#hmf0.update(n=0.966, sigma_8=0.829,cosmo_params={'Om0':0.315,'H0':67.74},Mmin=8,Mmax=11)
		nzbin = 81
		lz_base=np.linspace(0,40,nzbin)#np.array([ZT(np.log10(x)) for x in lt_base])
		#ln_M_z0 = np.zeros(nzbin)
		ln_M_z00 = np.zeros(nzbin)
		h = 0.6774
		for i in range(len(lz_base)):
			#hmf0.update(z=lz_base[i])
			#hmf0.update(delta_c = d_delta(lz_base[i]))#
			#lm = np.log10(hmf0.m/h)
			#ln = np.log10(hmf0.ngtm)#*h**3)
			#nm = interp1d(lm,ln)
			#ln_M_z0[i] = (10**nm(9)-10**nm(Mmax))

			hmf00.update(z=lz_base[i])
			hmf00.update(delta_c = d_delta(lz_base[i]))#
			lm = np.log10(hmf00.m/h)
			ln = np.log10(hmf00.ngtm)#*h**3)
			nm = interp1d(lm,ln)
			ln_M_z00[i] = (10**nm(-2*np.log10((1+lz_base[i])/10)+6)-10**nm(np.log10(2.5*((1+lz_base[i])/10)**-1.5)+7))

		hmf00.update(z=6)
		lm = np.log10(hmf00.m/h)
		ln = np.log10(hmf00.ngtm)#*h**3)
		nm = interp1d(lm,ln)
		ln_M_z_norm = (10**nm(-2*np.log10((1+6)/10)+6)-10**nm(np.log10(2.5*((1+6)/10)**-1.5)+7))
		totxt('Nps_z.txt',[lz_base, ln_M_z00],0,0,0)
		totxt('Nps_norm.txt',[[ln_M_z_norm]],0,0,0)
	else:
		data = np.array(retxt('Nps_z.txt',2,0,0))
		lz_base, ln_M_z00 = data[0], data[1]
		ln_M_z_norm = retxt('Nps_norm.txt',1,0,0)[0][0]

	lJ21_z = [10.5*((1+lz_base[i])/7)**-1.5*ln_M_z00[i]/ln_M_z_norm for i in range(len(lz_base))]
	J21_z = interp1d(lz_base,np.log10(lJ21_z))

	#nz0 = interp1d(lz_base,np.log10(ln_M_z0))
	nz00 = interp1d(lz_base,np.log10(ln_M_z00))
	#def n_a_CDM(a,h=0.6774):
	#	z = 1/a-1
	#	return 10**nz0(z)
	
	def n_a_CDM0(a,h=0.6774):
		z = 1/a-1
		return 10**nz00(z)

	#test = Jnu_final(7.5, 310, L = L_nu0)
	#test0 = Jnu_cosmic(7.5,L=L_nu0,nu=310,n_a=n_a_CDM,mode=1)
	#print(test)
	
	JHII_z = lambda z: 0.3*(z<=6) + 0.3*(z>6)*n_a_CDM0(1/(1+z))/n_a_CDM0(1/7)

	if tag==0:
		lz = np.linspace(19.9,1/(1-5e-2)-1,nbin)
		#"""
		lJ0 = [Jnu_final(z,L=L_nu0,nu=310,dndm=dndm0) for z in lz]
		lJ1 = [Jnu_final(z,L=L_nu1,nu=310,dndm=dndm1) for z in lz]
		totxt(rep0+'Jnuz.txt',[lz,lJ0,lJ1],0,0,0)
		#"""
		lnu = 10**np.linspace(np.log10(50),3,nbin)
		lJnu0 = [Jnu_final(zend,L=L_nu0,nu=x,dndm=dndm0) for x in lnu]
		lJnu1 = [Jnu_final(zend,L=L_nu1,nu=x,dndm=dndm1) for x in lnu]
		totxt(rep0+'Jnu.txt',[lnu,lJnu0,lJnu1],0,0,0)
		T0 = Tnu(310,Jnu_final(zend,L=L_nu0,nu=310,dndm=dndm0))
		T1 = Tnu(310,Jnu_final(zend,L=L_nu1,nu=310,dndm=dndm1))
		totxt(rep0+'T310.txt',[[T0, T1]],0,0,0)

		lnu_ = 10**np.linspace(3,np.log10(10e4),nbin)
		lJnu0_ = [Jnu_final(zend,L=L_nu0,nu=x,dndm=dndm0) for x in lnu_]
		lJnu1_ = [Jnu_final(zend,L=L_nu1,nu=x,dndm=dndm1) for x in lnu_]
		totxt(rep0+'Jnu_.txt',[lnu_,lJnu0_,lJnu1_],0,0,0)

	else:
		dataz = np.array(retxt(rep0+'Jnuz.txt',3,0,0))
		datanu = np.array(retxt(rep0+'Jnu.txt',3,0,0))
		datanu_ = np.array(retxt(rep0+'Jnu_.txt',3,0,0))
		dT = retxt(rep0+'T310.txt',1,0,0)[0]
		T0, T1 = dT[0], dT[1]
		lz, lJ0, lJ1 = dataz[0], dataz[1], dataz[2]
		lnu, lJnu0, lJnu1 = datanu[0], datanu[1], datanu[2]
		lnu_, lJnu0_, lJnu1_ = datanu_[0], datanu_[1], datanu_[2]
		

	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax1.plot(lnu, Tnu(lnu,np.array(lJnu1)), label=r'Structure formation, '+lmodel_[1],lw=1)
	ax1.plot(lnu, Tnu(lnu,np.array(lJnu0)), label=r'Structure formation, '+lmodel_[0],ls='--',lw=1)
	ax1.plot(lnu, Tnu(lnu,JHII_z(6)),'-.',label=r'Minihalo',lw=1)
	ax1.plot(lnu[lnu>0], Tnu(lnu[lnu>0],10**J21_z(1420/lnu[lnu>0]-1)), ls=':',color='r', label=r'21 cm emission',lw=1)
	lTnu_IGM = [Tnu(x,Jnu_bg(x)) for x in lnu]
	ax1.plot(lnu, lTnu_IGM, ls='-.',lw=2,color='g',label=r'ionized diffuse IGM')
	ax1.plot(lnu,Tnu_SKA(lnu),'k--',label=r'SKA',lw=2)
	ax1.fill_between(lnu,1e3*Tnu_sky_ff(lnu,-1),1e3*Tnu_sky_ff(lnu,1),facecolor='gray',label=r'$T_{\mathrm{ff}}^{\mathrm{G}}$',alpha=0.5)
	yup = np.max([160.0,np.max(Tnu(lnu,np.array(lJnu1))),np.max(Tnu(lnu,np.array(lJnu0)))])*1.05
	ax1.plot([1420/7,1420/7],[1e-4,yup],lw=0.5,color='k')
	ax1.fill_between([1420/7,1400],[1e-4,1e-4],[yup,yup],facecolor='gray',alpha=0.2)
	ax1.set_xlim(50,1000)
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
	print('Tnu at 310 MHz: {} (CDM) mK'.format(T0))
	print('Tnu at 310 MHz: {} (WDM) mK'.format(T1))

	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax1.plot(lnu_/1e3, Tnu(lnu_,np.array(lJnu1_)), label=r'Structure formation, '+lmodel_[1],lw=1)
	ax1.plot(lnu_/1e3, Tnu(lnu_,np.array(lJnu0_)), label=r'Structure formation, '+lmodel_[0],ls='--',lw=1)
	ax1.plot(lnu_/1e3, Tnu(lnu_,JHII_z(6)),'-.',label=r'Minihalo',lw=1)
	ax1.plot(lnu_/1e3, [Tnu(x,Jnu_bg(x)) for x in lnu_], ls='-.',lw=2,color='g',label=r'ionized diffuse IGM')
	ax1.plot(lnu_/1e3,Tnu_SKA(lnu_),'k--',label=r'SKA', lw=2)
	ax1.fill_between(lnu_/1e3,1e3*Tnu_sky_ff(lnu_,-1),1e3*Tnu_sky_ff(lnu_,1),facecolor='gray',label=r'$T_{\mathrm{ff}}^{\mathrm{G}}$',alpha=0.5)
	ax1.set_xlim(1.0,90)
	yup = np.max([0.31,np.max(Tnu(lnu_,np.array(lJnu1_))),np.max(Tnu(lnu_,np.array(lJnu0_)))])*1.05
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
	for i in range(len(lz)):
		if lz[i]<5:
			if lJ1[i] < lJ1[i-1]:
				lJ1[i]=lJ1[i-1]
			if lJ0[i] < lJ0[i-1]:
				lJ0[i]=lJ0[i-1]
	ax1.plot(lz, Tnu(310,np.array(lJ1)), label=r'$\nu_{\mathrm{obs}}=310\ \mathrm{MHz}$, '+lmodel_[1],lw=1)
	ax1.plot(lz, Tnu(310,np.array(lJ0)), '--',label=r'$\nu_{\mathrm{obs}}=310\ \mathrm{MHz}$, '+lmodel_[0],lw=1)
	yup = np.max([60, np.max(Tnu(310,np.array(lJ0))), np.max(Tnu(310,np.array(lJ1)))])*1.05
	ax1.plot([6,6],[1e-9,yup],lw=0.5,color='k')
	ax1.fill_between([0,6],[1e-9,1e-9],[yup,yup],facecolor='gray',alpha=0.2)
	ax1.set_xlim(0,20)
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


	"""
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
	"""


	"""
	lLp0 = retxt(rep0+'Lp_syn_cdm.txt',3,0,0)
	Lp0 = interp1d(lLp0[2],np.log10(meanL()*np.array(lLp0[1])))

	lLp1 = retxt(rep0+'Lp_syn_wdm.txt',3,0,0)
	Lp1 = interp1d(lLp1[2],np.log10(meanL()*np.array(lLp1[1])))

	p_index = 2.198
	A_syn0 = 10**Lp0(p_index)#1.1202483164633895e+59
	A_syn1 = 10**Lp1(p_index)#1.8152757384235942e+58
	facB = 1.0
	facn = 3.4619e-5*0.6774**3*2 * 0.75/0.7523662791819455 /(meanL()*2)
	print('beta_n = {}, fac = {}'.format(facn, meanL()))
	#facn = 1.0
	L_nu_syn0 = lambda x: jnu_syn_(10**x, A_syn0, p_index)*facB**((p_index+1)/4)*facn
	L_nu_syn1 = lambda x: jnu_syn_(10**x, A_syn1, p_index)*facB**((p_index+1)/4)*facn

	redshift = lLp0[0][0]
	dL = DZ(redshift)*(1+redshift)
	print('Syn Flux: {} [erg s^-1 cm^-2 Hz^-1]'.format((1+redshift)*L_nu_syn0(np.log10(1.4e9*(1+redshift)))/dL**2/4/np.pi))

	lnu = 10**np.linspace(np.log10(50),np.log10(1400),100)
	lJnu_syn0 = [Jnu_cosmic(zend,L=L_nu_syn0,nu=x,n_a=n_a_CDM,mode=tmode) for x in lnu]
	lJnu_syn1 = [Jnu_cosmic(zend,L=L_nu_syn1,nu=x,n_a=n_a_WDM,mode=tmode) for x in lnu]
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
	"""

	"""
	plt.figure()
	plt.plot(lz_base,np.array(ln_M_z1)/np.array(ln_M_z1_),label=lmodel[1])
	plt.plot(lz_base,np.array(ln_M_z0)/np.array(ln_M_z0_),'--',label=lmodel[0])
	plt.xlabel(r'$z$')
	plt.ylabel(r'Ratio')
	plt.legend()
	plt.xlim(7.5, 20)
	plt.ylim(10,1e4)
	plt.yscale('log')
	plt.tight_layout()
	plt.savefig('nM_ratio_z.pdf')
	plt.show()
	
	plt.figure()
	plt.plot(lz_base,ln_M_z00,label='EdS')
	plt.plot(lz_base,ln_M_z000,'--',label='Naoz et al. 2006')
	plt.xlabel(r'$z$')
	plt.ylabel(r'$N_{\mathrm{ps}}\ [h^{3}\mathrm{Mpc^{-3}}]$')
	plt.yscale('log')
	plt.legend()
	plt.tight_layout()
	plt.savefig(rep0+'nps_M_z.pdf')
	plt.show()
	
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


	lz_syn = np.linspace(3,30,100)#np.logspace(-2,np.log10(30),100)
	lJnu_syn0 = [Jnu_cosmic(max(x,zend),L=L_nu_syn0,nu=1420/(1+x),n_a=n_a_CDM,mode=tmode) for x in lz_syn]
	lJnu_syn1 = [Jnu_cosmic(max(x,zend),L=L_nu_syn1,nu=1420/(1+x),n_a=n_a_WDM,mode=tmode) for x in lz_syn]
	plt.figure()
	plt.plot(lz_syn, 2.725*(1+lz_syn)+Tnu(1420/(1+lz_syn),np.array(lJnu_syn1))/1e3,'--',label='Structure formation, '+lmodel[0])
	plt.plot(lz_syn, 2.725*(1+lz_syn)+Tnu(1420/(1+lz_syn),np.array(lJnu_syn0))/1e3,label='Structure formation, '+lmodel[1])
	plt.plot(lz_syn, 2.725*(1+lz_syn)+(Tnu_sky(1420/(1+lz_syn))-2.725)*0.075,'-.',label='7.5% of ARCADE 2 excess')
	plt.plot(lz_syn, 2.725*(1+lz_syn),':',label='no excess')
	plt.ylabel(r'$T\ [\mathrm{K}]$')
	plt.xlabel(r'$z$')
	plt.legend()
	plt.xlim(3,30)
	if mode==0:
		plt.yscale('log')
	plt.tight_layout()
	if mode==0:
		plt.savefig(rep0+'logTz_syn.pdf')
	else:
		plt.savefig(rep0+'Tz_syn.pdf')
	"""

	

