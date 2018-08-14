from onezone1 import *
from process import *

def phase_diagram(sn = 50, rep = './', indm = 0, edge = [0.005, 99.999], base = 'snapshot', ext = '.hdf5', mode = 0):
	ds = yt.load(rep+base+'_'+str(sn).zfill(3)+ext)
	ad = ds.all_data()
	keys = ds.field_list
	tag = np.sum([x[0] == 'PartType3' for x in keys])
	if tag>0:
		num_sink = len(ad[('PartType3','Masses')])
		print('Number of sink particles: {}'.format(num_sink))
	nd = ad[('PartType0','density')].to_equivalent("cm**-3", "number_density",mu=mmw(ad[('PartType0','Primordial HII')]))
	lT = temp(ad[('PartType0','InternalEnergy')],ad[('PartType0','Primordial HII')])
	lxh = ad[('PartType0','Primordial H2')]
	#lxhd = ad[('PartType0','Primordial HD')]
	#lxe = ad[('PartType0','Primordial e-')]
	#lMBE = MBE(lT, nd, lxe)
	#Mres = 32*ad[('PartType0','Masses')][0].to('Msun')

	#racc = np.array((Mres*UM/1e10/(4*np.pi*1e2*1.22*PROTON/3))**(1/3)/UL)
	#print('Accretion radius: {} kpc'.format(racc))

	rT = [0.99,6.5]#np.percentile(np.log10(lT), edge)
	rn = [-5.0,4.0]#np.percentile(np.log10(nd), edge)
	rx = np.percentile(np.log10(lxh), edge)
	#rxd = np.percentile(np.log10(lxhd), edge)
	#rxd[0] = max(rxd[0], -11)
	#rxe = np.percentile(np.log10(lxe), edge)
	#rMBE = np.percentile(np.log10(lMBE), edge)
	#print(np.max(nd))

	if mode==0:
		ref = main1(z=7.69, nini = n0(7.69), Xh2=1e-16)
		pT = np.log10(ref['T'])
		pn = np.log10(ref['n'])
		px = np.log10(ref['X'][3])
		totxt('old.txt',[pT,pn,px],0,0,0)
		ref0 = main1(Tini=6.7e4, z=7.69, nini = n0(7.69), Xh2=1e-16, mode=0)
		pT0 = np.log10(ref0['T'])
		pn0 = np.log10(ref0['n'])
		px0 = np.log10(ref0['X'][3])
		totxt('old0.txt',[pT0,pn0,px0],0,0,0)
	else:
		load = np.array(retxt('old.txt',3,0,0))
		pT = load[0]
		pn = load[1]
		px = load[2]
		load0 = np.array(retxt('old0.txt',3,0,0))
		pT0 = load0[0]
		pn0 = load0[1]
		px0 = load0[2]

	plt.figure()
	plt.subplot(111)
	plt.hist2d(np.log10(lT),np.log10(lxh),bins=100,norm=LogNorm(),range=[rT,rx])
	cb = plt.colorbar()
	cb.set_label(r'$\log(N)$')
	within = (rT[0]<pT)*(pT<rT[1])
	plt.plot(pT[within], px[within], 'r-',label='Minihalo')
	within0 = (rT[0]<pT0)*(pT0<rT[1])
	plt.plot(pT0[within0], px0[within0], 'r--',label='Shocked')
	plt.plot(np.log10([2.73*8.69,2.73*8.69]),rx,'k:',label='CMB')
	plt.legend()
	plt.ylabel(r'$\log([\mathrm{H_{2}/H}])$')
	plt.xlabel(r'$\log(T\ [\mathrm{K}])$')
	#plt.title(r'$[\mathrm{H_{2}/H}]-T$ phase diagram for '+lmodel[indm]+' at $z=$'+str(int(ds['Redshift']*100)/100),size=12)
	plt.tight_layout()
	plt.savefig(rep+'XH2T_'+lmodel[indm]+'_'+str(sn)+'.pdf')

	plt.figure()
	plt.subplot(111)
	plt.hist2d(np.log10(nd),np.log10(lxh),bins=100,norm=LogNorm(),range=[rn,rx])
	cb = plt.colorbar()
	cb.set_label(r'$\log(N)$')
	within = (rn[0]<pn)*(pn<rn[1])
	plt.plot(pn[within], px[within], 'r-',label='Minihalo')
	within0 = (rn[0]<pn0)*(pn0<rn[1])
	plt.plot(pn0[within0], px0[within0], 'r--',label='Shocked')
	plt.legend()
	plt.ylabel(r'$\log([\mathrm{H_{2}/H}])$')
	plt.xlabel(r'$\log(n\ [\mathrm{cm^{-3}}])$')
	#plt.title(r'$[\mathrm{H_{2}/H}]-n$ phase diagram for '+lmodel[indm]+' at $z=$'+str(int(ds['Redshift']*100)/100),size=12)
	plt.tight_layout()
	plt.savefig(rep+'XH2n_'+lmodel[indm]+'_'+str(sn)+'.pdf')

	plt.figure()
	plt.subplot(111)
	plt.hist2d(np.log10(nd),np.log10(lT),bins=100,norm=LogNorm(),range=[rn,rT])
	cb = plt.colorbar()
	cb.set_clim(1.0,1e6)
	cb.set_label(r'$\log(N)$')
	within = (rn[0]<pn)*(pn<rn[1])
	plt.plot(pn[within], pT[within], 'r-',label='Minihalo')
	within0 = (rn[0]<pn0)*(pn0<rn[1])
	plt.plot(pn0[within0], pT0[within0], 'r--',label='Shocked')
	plt.plot(rn,np.log10([2.73*8.69,2.73*8.69]),':',label='CMB',color='k')
	plt.legend()
	plt.xlabel(r'$\log(n\ [\mathrm{cm^{-3}}])$')
	plt.ylabel(r'$\log(T\ [\mathrm{K}])$')
	#plt.title(r'$T-n$ phase diagram for '+lmodel[indm]+' at $z=$'+str(int(ds['Redshift']*100)/100),size=12)
	plt.tight_layout()
	plt.savefig(rep+'Tn_'+lmodel[indm]+'_'+str(sn)+'.pdf')
	#plt.show()
	
	return ds

if __name__ == "__main__":
	sn = int(sys.argv[1])
	ds = phase_diagram(sn, indm=int(sys.argv[2]),mode=1)




