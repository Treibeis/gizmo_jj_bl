from txt import *
import yt
import caesar
from matplotlib.colors import LogNorm
from yt import YTQuantity
import matplotlib.patches as pltp
from matplotlib.collections import PatchCollection
from scipy.optimize import *
import sys
from cosmology import *

def temp_(u, Y, X = H_FRAC, gamma = 5.0/3):
	M = PROTON*4.0/(3*X+1+4*X*Y)
	U = u.to('erg/g')
	tem = M*U*(gamma-1)/BOL
	return tem

def temp(u, Y, X = H_FRAC, gamma = 5.0/3):
	M = PROTON*4.0/(3*X+1+4*X*Y)
	U = u*UE/UM
	tem = M*U*(gamma-1)/BOL
	return tem

def magneticf(u, rho, fac):
	U = u*UE/UM * rho
	B = (fac*U*8*np.pi)**0.5
	return B

def MBE(T, n, Y, X = H_FRAC, gamma = 5.0/3):
	mu = 4.0/(3*X+1+4*X*Y)
	xh = 4.0*X/(3*X+1)
	out = 1050*(T/200)**1.5/(mu/1.22)**2.0/(xh*n/1e4)**0.5 * (gamma/(5.0/3))**2
	return out
	
def mmw(xeH = XeH, xeHe = XeHe, X = Hfrac):
	xh = 4*X/(1+3*X)
	return 4.0/(1.0+3*X)/(xh*(1+xeH)+(1-xh)*(1+xeHe))

# dd["density"].to_equivalent("cm**-3", "number_density",mu=mmw())

def cosmicT(z = 99, zd = 200.0, zre = 1100):
	Tre = 2.73*(1+zre)
	if z>=zd:
		T = Tre*(1+z)/(1+zre)
	else:
		Td = Tre*(1+zd)/(1+zre)
		T = Td*(1+z)**2/(1+zd)**2
	return T
	

lls = ['-', '--', '-.', ':']
lls = lls*2
lrep = ['N128L4_cdm/', 'N128L4_wdm/', 'N128L4_cdm_dvx30/']
lmodel = ['CDM', 'WDM_3_kev', 'DVX30']

def loaddata(sn = 37, en = 50, rep = './', base = 'snapshot', ext = '.hdf5', post = 'caesar'):
	for i in range(sn, en+1):
		ds = yt.load(rep+base+'_'+str(i).zfill(3)+ext)
		obj = caesar.CAESAR(ds)
		obj.member_search()
		obj.save(rep+post+'_'+str(i).zfill(3)+ext)
	return {'obj':obj, 'ds':ds}

def recaesar(sn, rep = './', base = 'snapshot', ext = '.hdf5', post = 'caesar'):
	#ds = yt.load(rep+base+'_'+str(sn).zfill(3)+ext)
	obj = caesar.load(rep+post+'_'+str(sn).zfill(3)+ext)
	#obj.yt_dataset = ds
	return obj

def periodic(x, boxs = 4000):
	return x +  boxs * (-1*(x>boxs)+(x<0))

def virial_mass(group):
	return float(1e10*group.virial_quantities['circular_velocity']**2*group.radii['virial'].to('kpc/h')/G)

def to_center(pos0, ref, boxs = 4000):
	pos = np.copy(pos0)
	disp = boxs/2.0 - ref
	return periodic(pos+disp, boxs)

def center(pos0, ref, rad, boxs, fac = 2.5):
	pos = pos0 + boxs*(ref-pos0>boxs/2.0-rad*fac)*(ref>boxs/2.0) - boxs*(pos0-ref>boxs/2.0-rad*fac)*(ref<=boxs/2.0)
	return np.average(pos,axis=0)
	
def find_cen(group, boxs = 4000, h = 0.6774, fac = 2.5):
	rad = float(group.radii['virial'].to('kpccm/h'))
	DM = group.obj.data_manager
	ref = np.array(group.pos.to('kpccm/h')) 
	pos0 = DM.pos[DM.dmlist[group.dmlist]]*h
	return center(pos0, ref, rad, boxs, fac)#pos0 + boxs*(ref-pos0>boxs/2.0-rad*fac)*(ref>boxs/2.0) - boxs*(pos0-ref>boxs/2.0-rad*fac)*(ref<=boxs/2.0)
	#return np.average(pos,axis=0)

new_3Dax = lambda x: plt.subplot(x,projection='3d')

def dis3D(data, group_index, axis, nump = 1e5, norm = 10.0, boxs = 4000, h = 0.6774):
	group = data['obj'].halos[group_index]
	#ad = data['ds'].all_data()
	llb = [r'$x\ [h^{-1}\mathrm{kpc}]$', r'$y\ [h^{-1}\mathrm{kpc}]$', r'$z\ [h^{-1}\mathrm{kpc}]$']
	DM = group.obj.data_manager
	lm = virial_mass(group)
	rV = float(group.radii['virial'].to('kpccm/h'))
	posh = find_cen(group, boxs, h)
	posg = DM.pos[DM.glist[group.glist]]*h#.to('kpccm/h')
	posd = DM.pos[DM.dmlist[group.dmlist]]*h#.to('kpccm/h')
	posg = to_center(posg,posh,boxs).T
	posd = to_center(posd,posh,boxs).T
	lpdm = np.random.choice(posd.shape[1], min(int(nump),posd.shape[1]),replace=False)
	lpg = np.random.choice(posg.shape[1], min(int(nump),posg.shape[1]),replace=False)
	#axis = plt.subplot(111,projection='3d')
	#axis.set_aspect('equal','box')
	axis.plot(posd[0][lpdm],posd[1][lpdm],posd[2][lpdm],'.',markersize=0.1,label='Dark matter')
	axis.plot(posg[0][lpg],posg[1][lpg],posg[2][lpg],'.',markersize=0.1,label='Gas')
	axis.scatter(boxs/2,boxs/2,boxs/2,s=rV*norm,edgecolors='g',facecolors='g',lw=2,alpha=0.5)
	axis.set_xlabel(llb[0])
	axis.set_ylabel(llb[1])
	axis.set_zlabel(llb[2])
	return [posh, rV, lm]
	
	#print(np.max(DM.pos),np.min(DM.pos))
"""
def candidate(data, nh = 5, indm = 1, icf = 'ics_N128L4_wdm.dat', icrep = 'ics/', boxs = 4000):
	ic_ds = yt.load(icrep+icf)
	for i in range(nh):
		ax = new_3Dax(111)
		out = dis3D(data, i, ax)
		print('Mass centre position: {}, virial radius: {}({}), M [10^9 Msun]: {}'.format(out[0]/boxs, out[1], out[1]/boxs, out[2]/1e9))
		plt.show()
		tag = input('Do you want this halo? (y/n)')
		if tag=='y':
			data['obj'].halos[i].write_IC_mask(ic_ds, 'mask_'+str(i)+'_'+lmodel[indm])
"""

def select(data, indm = 1, nh = 3, icf = 'ics_N128L4_wdm.dat', icrep = 'ics/', M_rg = [1e9, 1e10], in_box = 1.0, Dres = 4, boxs = 4000, h = 0.6774):
	ic_ds = yt.load(icrep+icf)
	lh = data['obj'].halos
	h_info = []
	count = 0
	for i in range(len(lh)):
		MV = virial_mass(lh[i])
		#print(MV/1e9)
		if MV>M_rg[0] and MV<M_rg[1] and count<nh:
			posh = find_cen(lh[i], boxs, h)/boxs
			tag = [(posh[i]<=0.5+in_box)and(posh[i]>0.5-in_box) for i in range(3)]
			tag = tag[0] and tag[1] and tag[2]
			if tag==True:
				#data['obj'].halos[i].write_IC_mask(ic_ds, 'mask_'+str(i)+'_'+lmodel[indm])
				ax = new_3Dax(111)
				out = dis3D(data, i, ax)
				print('Mass center position: {}, virial radius: {}({}), M [10^9 Msun]: {}'.format(out[0]/boxs, out[1], out[1]/boxs, out[2]/1e9))
				plt.savefig('dis'+'_'+str(i)+'_'+lmodel[indm]+'.png')
				plt.show()
				Rtb = 200**(1/3)*(1.5*Dres+1)*out[1]/boxs
				print('Rtb: {}'.format(Rtb))
				#h_info.append([out[0][i]/boxs for i in range(3)]+[Rtb, out[2]/1e9])
				tag = input('Do you want this halo? (y/n)')
				if tag=='y':
					count += 1
					data['obj'].halos[i].write_IC_mask(ic_ds, 'mask_'+str(i)+'_'+lmodel[indm])
					mask = retxt('mask_'+str(i)+'_'+lmodel[indm], 3, 0, 0)
					#ax = new_3Dax(111)
					#ax.plot(mask[0],mask[1],mask[2],'.',markersize=0.1)
					#plt.show()
					mask = np.array(mask).T
					new_posh = center(mask, posh, out[1]/boxs, 1.0)
					mask = to_center(mask, new_posh, 1.0)
					ext = np.max(mask,axis=0)-np.min(mask, axis=0)
					reg_center = (np.max(mask,axis=0)+np.min(mask, axis=0))/2.0 + new_posh - 0.5
					h_info.append([i]+[reg_center[i] for i in range(3)]+[Rtb, out[2]/1e9]+[ext[i] for i in range(3)])
					print('reg_center: {}, reg_extent: {}'.format(new_posh, ext))
	totxt('halo_info_'+lmodel[indm],h_info,0,0,0)
	
import hmf
import hmf.wdm

hmf0 = hmf.MassFunction()
hmf0.update(n=0.966, sigma_8=0.829,cosmo_params={'Om0':0.315,'H0':67.74},Mmin=4,Mmax=11)
hmf1 = hmf.wdm.MassFunctionWDM()
hmf1.update(n=0.966, sigma_8=0.829,cosmo_params={'Om0':0.315,'H0':67.74},Mmin=4,Mmax=11,wdm_mass=3)

#HMF = 'mVector_PLANCK-z8.8 .txt'
#HMF = 'mVector_PLANCK-SMT .txt'
#HMF0 = 'mVector_PLANCK-SMT0 .txt'
Redf=retxt('output_25_4.txt',1,0,0)[0]
fac0 = 0.208671*0.3153754*0.1948407



def halom(sn, rep ='./', indm = 0, anaf1 = hmf1, anaf0 = hmf0, mlow = 1e6, fac = fac0, edge = [0.0, 99.5], nbin = 15, boxs = 4.0 ,base = 'snapshot', ext = '.hdf5', post = 'caesar'):
	#ds = yt.load(rep+base+'_'+str(sn).zfill(3)+ext)
	z = 1/Redf[sn]-1
	obj = caesar.load(rep+post+'_'+str(sn).zfill(3)+ext)
	lh = obj.halos
	#a = 1/(ds['Redshift']+1)
	lm0 = np.array([x.masses['total'] for x in lh if x.masses['total']>=mlow])
	lm = np.array([(1e10*x.virial_quantities['circular_velocity']**2*x.radii['virial'].to('kpc/h')/G) for x in lh if x.masses['total']>=mlow])
	rm0 = np.percentile(np.log10(lm0), edge)
	rm = np.percentile(np.log10(lm), edge)
	his1, ed1 = np.histogram(np.log10(lm0), nbin, rm0)
	b1 = (ed1[1:]+ed1[:-1])/2
	bs1 = ed1[1]-ed1[0]
	his2, ed2 = np.histogram(np.log10(lm), nbin, rm)
	b2 = (ed2[1:]+ed2[:-1])/2
	bs2 = ed2[1]-ed2[0]
	anaf1.update(z=z)
	lm1 = np.log10(anaf1.m)
	ln1 = np.log10(anaf1.dndlog10m)
	mf_ana1 = interp1d(lm1, ln1)
	anaf0.update(z=z)
	lm0 = np.log10(anaf0.m)
	ln0 = np.log10(anaf0.dndlog10m)
	mf_ana0 = interp1d(lm0, ln0)
	plt.figure()
	plt.errorbar(b1, his1/boxs**3/bs1/fac, yerr=his1**0.5/boxs**3/bs1/fac, label='total mass')
	plt.errorbar(b2, his2/boxs**3/bs2/fac, yerr=his2**0.5/boxs**3/bs2/fac, ls='--', label='virial mass')
	plt.plot(b2, 10**mf_ana1(b2),'-.', label='Sheth-Mo-Tormen (WDM)')
	plt.plot(b2, 10**mf_ana0(b2),':', label='Sheth-Mo-Tormen (CDM)')
	plt.legend()
	plt.yscale('log')
	plt.title(r'Halo mass function for '+lmodel[indm]+' at $z=$'+str(int(z*100)/100),size=14)
	plt.xlabel(r'$\log(M\ [h^{-1}M_{\odot}])$')
	plt.ylabel(r'$\frac{dn}{d\log(M)}\ [h^{3}\mathrm{Mpc^{-3}}]$')
	plt.tight_layout()
	plt.savefig(rep+'hmass_'+lmodel[indm]+'_'+str(sn)+'.pdf')
	#plt.show()
	return lm

def plothm(l = lrep, sn = 28, nbin = 15, rm = [8, 10], anaf = hmf0, boxs = 4.0):
	z = 1/Redf[sn]-1
	lhm = [halom(sn, l[x], x) for x in range(3)]
	lhis = []
	lbs = []
	for i in range(3):
		his, ed = np.histogram(np.log10(lhm[i]), nbin, rm)
		lbs.append((ed[1:]+ed[:-1])/2+0.03*i)
		if i==0: 
			bs = (ed[1]-ed[0])
		lhis.append(his)
	anaf.update(z=z)
	lm = np.log10(anaf.m)
	ln = np.log10(anaf.dndlog10m)
	mf_ana = interp1d(lm, ln)
	#mf_file = retxt(anaf, 11, 12, 0)
	#mf_ana = interp1d(np.log10(mf_file[0]),mf_file[7])
	plt.figure(figsize=(12,5))
	plt.subplot(121)
	[plt.errorbar(lbs[i], lhis[i]/boxs**3/bs, yerr=lhis[i]**0.5/boxs**3/bs, ls=lls[i], label=lmodel[i]) for i in range(3)]
	plt.plot(lbs[0], mf_ana(lbs[0]), ls=lls[3], label='Sheth-Mo-Tormen')
	plt.legend()
	plt.yscale('log')
	plt.xlabel(r'$\log(M_{V}\ [h^{-1}M_{\odot}])$')
	plt.ylabel(r'$\frac{dn}{d\log(M)}\ [h^{3}\mathrm{Mpc^{-3}}]$')
	plt.subplot(122)
	plt.errorbar(lbs[0], mf_ana(lbs[0])-lhis[0]/boxs**3/bs, yerr=lhis[0]**0.5/boxs**3/bs, ls=lls[0], label='Analytical$-$CDM')
	[plt.errorbar(lbs[i],(lhis[i]-lhis[0])/boxs**3/bs, yerr=(lhis[i]+lhis[0])**0.5/boxs**3/bs, ls =lls[i], label=lmodel[i]+r'$-$'+lmodel[0]) for i in range(1,3)]
	plt.plot(lbs[0],np.zeros(lbs[0].shape),color='gray',lw=2.0,alpha=0.5)
	plt.legend()
	plt.xlabel(r'$\log(M_{V}\ [h^{-1}M_{\odot}])$')
	plt.ylabel(r'$\Delta\frac{dn}{d\log(M)}\ [h^{3}\mathrm{Mpc^{-3}}]$')
	plt.tight_layout()
	#plt.suptitle(r'Halo mass function at $z=$'+str(int(z*100)/100),size=14)
	plt.savefig('dhmass_'+str(sn)+'.pdf')
	plt.show()
		
rec_default = [[1,1,1]]*1+[np.array([0.3783187, 0.3497985, 0.399293243])*4000]

#[[500,400],[900,750],[500,200],[500,400]]+[[150,340],[300,500],[1,1],[500,400],[300,600],[1,1],[1,1],[700,550]]#[[200,200] for i in range(10)]	

def cosweb(sn, indm, data, rep = './', ax = [0, 1], rec0 = rec_default, mode = 0, norm = 1e6, nump = 1e5, nb = 100, boxs = 4000, h = 0.6774, norm2 = 50):
	llb = [r'$x\ [h^{-1}\mathrm{kpc}]$', r'$y\ [h^{-1}\mathrm{kpc}]$', r'$z\ [h^{-1}\mathrm{kpc}]$']
	#data = loaddata(sn, sn, rep, base, ext, post)
	ds = data['ds']#yt.load(rep+base+'_'+str(sn).zfill(3)+ext)
	ad = ds.all_data()
	obj = data['obj']#caesar.load(rep+post+'_'+str(sn).zfill(3)+ext)
	lh = obj.halos
	posh = np.array([find_cen(x, boxs, h) for x in lh]).T
	#lm = np.array([x.masses['total'] for x in lh])
	lm = np.array([(1e10*x.virial_quantities['circular_velocity']**2*x.radii['virial'].to('kpc/h')/G) for x in lh])
	posd = np.array(ad[('PartType1','Coordinates')].to('kpccm/h')).T
	posg = np.array(ad[('PartType0','Coordinates')].to('kpccm/h')).T
	if mode==-1:
		poss = np.array(ad[('PartType3','Coordinates')].to('kpccm/h')).T
		lms = np.array(ad[('PartType3','Masses')]/ad[('PartType0','Masses')][0])/norm2
		#print(lms)
	if rec0!=[]:
		rec = np.array(rec0)
		lrp = [find_cen(lh[x], boxs, h) for x in range(len(rec0))]
		#print(lrp)
	print('Mass resolution:',ad[('PartType0','Masses')][0].to('Msun/h'), ad[('PartType1','Masses')][0].to('Msun/h'))
	print('Binsize:',ds['BoxSize']/nb,'kpccm/h')
	print('Number of DM halos:',posh.shape[1])
	print('Halo mass list {}'.format([np.log10(lm[i]) for i in range(len(rec0))]))
	if mode==1:
		lpdm = np.random.choice(posd.shape[1], min(int(nump),posd.shape[1]),replace=False)
		lpg = np.random.choice(posg.shape[1], min(int(nump),posg.shape[1]),replace=False)
		axis = plt.subplot(111,projection='3d')
		axis.plot(posd[0][lpdm],posd[1][lpdm],posd[2][lpdm],'.',markersize=0.1,label='Dark matter')
		axis.plot(posg[0][lpg],posg[1][lpg],posg[2][lpg],'.',markersize=0.1,label='Gas')
		axis.scatter(posh[0],posh[1],posh[2],s=lm/norm,edgecolors='g',facecolors='g',lw=2,alpha=0.5)
		axis.set_xlabel(llb[0])
		axis.set_ylabel(llb[1])
		axis.set_zlabel(llb[2])
		#axis.set_xlim(0,4000)
		#axis.set_ylim(0,4000)
		#axis.set_zlim(0,4000)
		axis.set_title(r'Particle distribution for '+lmodel[indm]+' at $z=$'+str(int(ds['Redshift']*100)/100),size=14)
		plt.tight_layout()
		plt.savefig(rep+'DM_dis_'+lmodel[indm]+'_'+str(posh.shape[1])+'_'+str(sn)+'.png')
		#plt.show()

	plt.figure(figsize=(12,5))
	ax1 = plt.subplot(121)
	plt.hist2d(posd[ax[0]],posd[ax[1]],bins=nb,norm=LogNorm(),cmap=plt.cm.coolwarm)
	cb = plt.colorbar()
#	cb.set_label(r'$\log(N)$')
	if rec0!=[]:
		lrec = [pltp.Rectangle((lrp[i]-rec[i]/2)[ax],rec[i][ax[0]],rec[i][ax[1]],fill=False) for i in range(len(rec0))]
		pc = PatchCollection(lrec,alpha=0.5)
		ax1.add_collection(pc)
	#plt.scatter(posh[ax[0]],posh[ax[1]],s=lm/norm,edgecolors='g',facecolors='none',lw=1)
	if mode==-1:
		plt.scatter(poss[ax[0]],poss[ax[1]], s=lms, edgecolors='purple', facecolors='none',lw=1)
	plt.xlabel(llb[ax[0]])
	plt.ylabel(llb[ax[1]])
	plt.title(r'Dark matter',size=14)
#	plt.tight_layout()
	ax2 = plt.subplot(122)
	plt.hist2d(posg[ax[0]],posg[ax[1]],bins=nb,norm=LogNorm(),cmap=plt.cm.hot)
	cb = plt.colorbar()
	#cb.set_label(r'$\log(N)$')
	if rec0!=[]:
		lrec = [pltp.Rectangle((lrp[i]-rec[i]/2)[ax],rec[i][ax[0]],rec[i][ax[1]],fill=False) for i in range(len(rec0))]
		pc = PatchCollection(lrec,alpha=0.5,facecolor='white')
		ax2.add_collection(pc)
	#plt.scatter(posh[ax[0]],posh[ax[1]],s=lm/norm,edgecolors='g',facecolors='none',lw=1)
	if mode==-1:
		plt.scatter(poss[ax[0]],poss[ax[1]], s=lms, edgecolors='purple', facecolors='none',lw=1)
	plt.xlabel(llb[ax[0]])
	plt.ylabel(llb[ax[1]])
	plt.title(r'Gas',size=14)
	plt.suptitle(r'Cosmic web for '+lmodel[indm]+' at $z=$'+str(int(ds['Redshift']*100)/100),size=14)
	plt.tight_layout()
	plt.savefig(rep+'cswb_'+lmodel[indm]+'_'+str(ax[0])+str(ax[1])+'_'+str(sn)+'.pdf')
	plt.show()
	return ds
	

def phase(sn = 50, rep = './', indm = 0, edge = [0.01, 100.0], base = 'snapshot', ext = '.hdf5'):
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
	lxhd = ad[('PartType0','Primordial HD')]
	lxe = ad[('PartType0','Primordial e-')]
	lMBE = MBE(lT, nd, lxe)
	Mres = 32*ad[('PartType0','Masses')][0].to('Msun')

	racc = np.array((Mres*UM/1e10/(4*np.pi*1e2*1.22*PROTON/3))**(1/3)/UL)
	print('Accretion radius: {} kpc'.format(racc))

	rT = [0.99,6.5]#np.percentile(np.log10(lT), edge)
	rn = [-5.0,4.0]#np.percentile(np.log10(nd), edge)
	rx = np.percentile(np.log10(lxh), edge)
	rxd = np.percentile(np.log10(lxhd), edge)
	rxd[0] = max(rxd[0], -11)
	rxe = np.percentile(np.log10(lxe), edge)
	rMBE = np.percentile(np.log10(lMBE), edge)
	print(np.max(nd))

	plt.figure()
	plt.subplot(111)
	plt.hist2d(np.log10(lT),np.log10(lxh),bins=100,norm=LogNorm(),range=[rT,rx])
	cb = plt.colorbar()
	cb.set_label(r'$\log(N)$')
	plt.ylabel(r'$\log([\mathrm{H_{2}/H}])$')
	plt.xlabel(r'$\log(T\ [\mathrm{K}])$')
	plt.title(r'$[\mathrm{H_{2}/H}]-T$ phase diagram for '+lmodel[indm]+' at $z=$'+str(int(ds['Redshift']*100)/100),size=12)
	plt.tight_layout()
	plt.savefig(rep+'XH2T_'+lmodel[indm]+'_'+str(sn)+'.pdf')

	plt.figure()
	plt.subplot(111)
	plt.hist2d(np.log10(nd),np.log10(lT),bins=100,norm=LogNorm(),range=[rn,rT])
	plt.plot(rn,np.log10([2.73*ds['Redshift'],2.73*ds['Redshift']]),'k:',label='CMB')
	plt.legend()
	cb = plt.colorbar()
	cb.set_label(r'$\log(N)$')
	cb.set_clim(1.0,1e6)
	plt.xlabel(r'$\log(n\ [\mathrm{cm^{-3}}])$')
	plt.ylabel(r'$\log(T\ [\mathrm{K}])$')
	#plt.title(r'$T-n$ phase diagram for '+lmodel[indm]+' at $z=$'+str(int(ds['Redshift']*100)/100),size=12)
	plt.tight_layout()
	plt.savefig(rep+'Tn_'+lmodel[indm]+'_'+str(sn)+'.pdf')
	#plt.show()
	plt.figure()
	plt.subplot(111)
	plt.hist2d(np.log10(nd),np.log10(lxh),bins=100,norm=LogNorm(),range=[rn,rx])
	cb = plt.colorbar()
	cb.set_label(r'$\log(N)$')
	plt.ylabel(r'$\log([\mathrm{H_{2}/H}])$')
	plt.xlabel(r'$\log(n\ [\mathrm{cm^{-3}}])$')
	plt.title(r'$[\mathrm{H_{2}/H}]-n$ phase diagram for '+lmodel[indm]+' at $z=$'+str(int(ds['Redshift']*100)/100),size=12)
	plt.tight_layout()
	plt.savefig(rep+'XH2n_'+lmodel[indm]+'_'+str(sn)+'.pdf')
	#plt.show()
	plt.figure()
	plt.subplot(111)
	plt.hist2d(np.log10(nd),np.log10(lxhd),bins=100,norm=LogNorm(),range=[rn,rxd])
	cb = plt.colorbar()
	cb.set_label(r'$\log(N)$')
	plt.ylabel(r'$\log([\mathrm{HD/D}])$')
	plt.xlabel(r'$\log(n\ [\mathrm{cm^{-3}}])$')
	plt.title(r'$[\mathrm{HD/D}]-n$ phase diagram for '+lmodel[indm]+' at $z=$'+str(int(ds['Redshift']*100)/100),size=12)
	plt.tight_layout()
	plt.savefig(rep+'XHDn_'+lmodel[indm]+'_'+str(sn)+'.pdf')
	#plt.show()
	plt.figure()
	plt.subplot(111)
	plt.hist2d(np.log10(nd),np.log10(lxe),bins=100,norm=LogNorm(),range=[rn,rxe])
	cb = plt.colorbar()
	cb.set_label(r'$\log(N)$')
	plt.ylabel(r'$\log([\mathrm{e^{-}/H}])$')
	plt.xlabel(r'$\log(n\ [\mathrm{cm^{-3}}])$')
	plt.title(r'$[\mathrm{e^{-}/H}]-n$ phase diagram for '+lmodel[indm]+' at $z=$'+str(int(ds['Redshift']*100)/100),size=12)
	plt.tight_layout()
	plt.savefig(rep+'Xen_'+lmodel[indm]+'_'+str(sn)+'.pdf')
	#plt.show()
	plt.figure()
	plt.subplot(111)
	plt.hist2d(np.log10(lT),np.log10(lxe),bins=100,norm=LogNorm(),range=[rT,rxe])
	cb = plt.colorbar()
	cb.set_label(r'$\log(N)$')
	plt.ylabel(r'$\log([\mathrm{e^{-}/H}])$')
	plt.xlabel(r'$\log(T\ [\mathrm{K}])$')
	plt.title(r'$[\mathrm{e^{-}/H}]-T$ phase diagram for '+lmodel[indm]+' at $z=$'+str(int(ds['Redshift']*100)/100),size=12)
	plt.tight_layout()
	plt.savefig(rep+'XeT_'+lmodel[indm]+'_'+str(sn)+'.pdf')
	
	plt.figure()
	plt.subplot(111)
	plt.hist2d(np.log10(nd),np.log10(lMBE),bins=100,norm=LogNorm(),range=[rn,rMBE])
	cb = plt.colorbar()
	cb.set_label(r'$\log(N)$')
	plt.plot(rn,np.log10([Mres,Mres]),label=r'$M_{\mathrm{res}}$',color='k')
	plt.ylabel(r'$\log(M_{\mathrm{BE}}\ [M_{\odot}])$')
	plt.xlabel(r'$\log(n\ [\mathrm{cm^{-3}}])$')
	plt.title(r'$M_{\mathrm{BE}}-n$ phase diagram for '+lmodel[indm]+' at $z=$'+str(int(ds['Redshift']*100)/100),size=12)
	plt.tight_layout()
	plt.savefig(rep+'MBEn_'+lmodel[indm]+'_'+str(sn)+'.pdf')
	#plt.show()
	
	return ds

def sinks(ntask = 4, rep = 'sink/', base = 'sink', ext = '.txt', mode=0):
	raw = []
	for i in range(ntask):
		d = retxt(rep+base+'.'+str(i)+ext,10,0,0)
		raw.append(d)
	raw = np.hstack(raw).T
	raw = np.array(sorted(raw, key=lambda x: x[0])).T
	if mode!=0:
		from radio import TZ
		raw[0] = np.array([TZ(1/x-1) for x in raw[0]])
		raw[1] = np.array([TZ(1/x-1) for x in raw[1]])
	nsink = int(np.max(raw[9]))+1
	out = []
	for i in range(nsink):
		out.append([x[raw[9]==i] for x in raw])
	return out

if __name__ == "__main__":
	sn = int(sys.argv[1])
	if len(sys.argv)<3:
		model = 1
	else:
		model = int(sys.argv[2])
	#data = loaddata(sn,sn)
	#ds = cosweb(sn, model, data, mode = -1, norm = 1e6, nb = 500)
	#ds = cosweb(sn, model, data, ax = [1,2], norm = 1e6, mode = -1, nb = 500)
	#ds = cosweb(sn, model, data, ax = [2,0], norm = 1e6, mode = -1, nb = 500)
	ds = phase(sn,indm=model)
	fac0 = 1.0 #0.3783187*0.3497985*0.399293243
	#lm = halom(sn, indm = model, fac = fac0)
	#plothm(sn=19)

	"""
	from radio import *
	lmst = np.array(retxt('sinkfr.txt',3,0,0))
	ntask = 160
	dsink = sinks(ntask,mode=1)
	lt = np.array([x[1][0]/1e9/YR for x in dsink])
	lm0 = np.array([x[2][0] for x in dsink])
	lm1 = np.array([x[2][-1] for x in dsink])
	
	plt.figure()
	plt.scatter(lt, lm0, label='initial',alpha=0.3)
	plt.scatter(lt, lm1, label='final',alpha=0.3)
	plt.yscale('log')
	plt.xlabel(r'$t_{0}\ [\mathrm{Gyr}]$')
	plt.ylabel(r'$m_{\mathrm{sink}}\ [M_{\odot}]$')
	plt.legend()
	plt.tight_layout()
	plt.savefig('msink_t0_'+str(lmodel[model])+'.pdf')

	plt.figure()
	plt.scatter(lt, lm1/lm0, label='initial',alpha=0.6)
	plt.yscale('log')
	plt.xlabel(r'$t_{0}\ [\mathrm{Gyr}]$')
	plt.ylabel(r'$m_{\mathrm{sink},1}/m_{\mathrm{sink},0}$')
	plt.tight_layout()
	plt.savefig('msink_ratio_t0_'+str(lmodel[model])+'.pdf')

	plt.figure()
	a = plt.hist(lt,density=True,bins=50,alpha=0.5)
	plt.xlabel(r'$t_{0}\ [\mathrm{Gyr}]$')
	plt.ylabel(r'Probability density')
	plt.tight_layout()
	plt.savefig('t0_'+str(lmodel[model])+'.pdf')

	plt.figure(figsize=(10,4))
	plt.subplot(121)
	plt.plot(np.array([TZ(1/x-1)/1e9/YR for x in lmst[0]])[lmst[1]>0], lmst[2][lmst[1]>0])
	plt.xlabel(r'$t_{U}\ [\mathrm{Gyr}]$')
	plt.ylabel(r'$M_{\mathrm{sink}}\ [M_{\odot}]$')
	#plt.tight_layout()
	#plt.savefig('Msink_tU'+str(lmodel[model])+'.pdf')
	plt.subplot(122)
	plt.plot((1/lmst[0]-1)[lmst[1]>0], lmst[2][lmst[1]>0])
	plt.xlabel(r'$z$')
	plt.ylabel(r'$M_{\mathrm{sink}}\ [M_{\odot}]$')
	plt.tight_layout()
	plt.savefig('Msink_tU_z_'+str(lmodel[model])+'.pdf')

	plt.figure(figsize=(10,4))
	plt.subplot(121)
	plt.plot(np.array([TZ(1/x-1)/1e9/YR for x in lmst[0]])[lmst[1]>0], lmst[1][lmst[1]>0])
	plt.xlabel(r'$t_{U}\ [\mathrm{Gyr}]$')
	plt.ylabel(r'$N_{\mathrm{sink}}$')
	#plt.tight_layout()
	#plt.savefig('Msink_tU'+str(lmodel[model])+'.pdf')
	plt.subplot(122)
	plt.plot((1/lmst[0]-1)[lmst[1]>0], lmst[1][lmst[1]>0])
	plt.xlabel(r'$z$')
	plt.ylabel(r'$N_{\mathrm{sink}}$')
	plt.tight_layout()
	plt.savefig('Nsink_tU_z_'+str(lmodel[model])+'.pdf')

	plt.show()
	"""

"""
def plot_Tn(sn = 25, en = 50, N = 4, base = 'snapshot', post = 'caesar', ext = '.hdf5'):
	X = Hfrac
	mu = PROTON*4.0/(1+3*X)
	lhalo = []
	lT = []
	ln = []
	lm = []
	lz = []
	for i in range(sn, en+1):
		ds = yt.load(base+'_'+str(i).zfill(3)+ext)
		obj = caesar.load(post+'_'+str(i).zfill(3)+ext)
		for h in obj.halos:
			data = [h.local_mass_density.to('g/cm**3')/mu, h.temperature, h.masses['total'], ds['Redshift']]
			lhalo.append(data)
			lT.append([obj.halos[i].temperature for i in range(N)])
			ln.append([obj.halos[i].local_mass_density.to('g/cm**3')/mu for i in range(N)])
			lm.append([obj.halos[i].masses['total'] for i in range(N)])
			lz.append(ds['Redshift'])
	lhalo = np.array(np.matrix(lhalo).transpose())
	lT = np.array(np.matrix(lT).transpose())
	ln = np.array(np.matrix(ln).transpose())
	lm = np.array(np.matrix(lm).transpose())

	plt.figure(figsize=(12,6))
	plt.subplot(121)
	plt.scatter(lhalo[0],lhalo[1],c=lhalo[3],cmap=plt.cm.coolwarm,s=20)
	cb = plt.colorbar()
	cb.ax.set_ylabel('Redshift')
	plt.xlim(min(lhalo[0])/3.0,max(lhalo[0])*3.0)
	plt.xlabel(r'$\bar{n}\ [\mathrm{cm^{-3}}]$')
	plt.ylabel(r'$\bar{T}\ [\mathrm{K}]$')
	plt.yscale('log')
	plt.xscale('log')
	plt.subplot(122)
	plt.scatter(lhalo[0],lhalo[1],c=np.log10(lhalo[2]),cmap=plt.cm.jet,s=20)
	cb = plt.colorbar()
	cb.ax.set_ylabel(r'$\log(M/M_{\odot})$')
	plt.xlim(min(lhalo[0])/3.0,max(lhalo[0])*3.0)
	plt.xlabel(r'$\bar{n}\ [\mathrm{cm^{-3}}]$')
	#plt.ylabel(r'$T\ [\mathrm{K}]$')
	plt.yscale('log')
	plt.xscale('log')
	plt.suptitle(r'Distribution of DM halos in $\bar{T}-\bar{n}$ space')
	#plt.tight_layout()
	plt.savefig('Tn1.pdf')
	plt.show()

	plt.figure(figsize=(12,6))
	plt.subplot(121)
	b = [plt.plot(ln[i], lT[i], ls = lls[i], label=r'Halo_'+str(i),color='gray',lw=0.8) for i in range(int(N/2))]
	plt.legend()
	a = [plt.scatter(ln[i], lT[i], c=lz, cmap=plt.cm.coolwarm, alpha=0.5) for i in range(int(N/2))]
	cb = plt.colorbar()
	cb.ax.set_ylabel('Redshift')
	plt.xlabel(r'$\bar{n}\ [\mathrm{cm^{-3}}]$')
	plt.ylabel(r'$\bar{T}\ [\mathrm{K}]$')
	plt.yscale('log')
	plt.xscale('log')
	plt.xlim(min([min(ln[i])/3.0 for i in range(int(N/2))]),max([max(ln[i])*3.0 for i in range(int(N/2))]))
	plt.subplot(122)
	b = [plt.plot(ln[i], lT[i], ls = lls[i], label=r'Halo_'+str(i),color='gray',lw=0.8) for i in range(int(N/2),N)]
	plt.legend()
	a = [plt.scatter(ln[i], lT[i], c=np.log10(lm[i]), cmap=plt.cm.jet, alpha=0.5) for i in range(int(N/2),N)]
	cb = plt.colorbar()
	cb.ax.set_ylabel(r'$\log(M/M_{\odot})$')
	plt.xlabel(r'$\bar{n}\ [\mathrm{cm^{-3}}]$')
	#plt.ylabel(r'$T\ [\mathrm{K}]$')
	plt.yscale('log')
	plt.xscale('log')
	#plt.ylim(100,1000)
	plt.xlim(min([min(ln[i])/3.0 for i in range(int(N/2),N)]),max([max(ln[i])*3.0 for i in range(int(N/2),N)]))
	plt.tight_layout()
	plt.suptitle('Evolution of the most massive 4 halos (NOT exactly tracked with the merger-tree)',size=12)
	plt.savefig('Tn2.pdf')
	plt.show()
"""		
		
