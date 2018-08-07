from process import *
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.misc import derivative
from sklearn.neighbors import NearestNeighbors
import time
import multiprocessing as mp
from matplotlib.colors import LogNorm
import scipy.stats as stats
from cosmology import *

def T_cosmic(z, alpha = -4, beta = 1.27, z0 = 189.6, zi = 1100, T0 = 3300):
	def integrand(logt):
		return alpha/3.0-(2+alpha)/3.0*(1-np.exp(-(ZT(logt)/z0)**beta))
	I = quad(integrand, np.log10(TZ(zi)/1e9/YR), np.log10(TZ(z)/1e9/YR))[0]#, epsrel = 1e-6)[0]
	temp = T0*10**I
	return temp

def rhom(a, Om = 0.315, h = 0.6774):
	H0 = h*100*UV/UL/1e3
	rho0 = Om*H0**2*3/8/np.pi/GRA
	return rho0/a**3

def tff(z = 10.0, delta = 200):
	return (3*np.pi/(32*GRA*delta*rhom(1/(1+z))))**0.5

def Lvir(m = 1e10, z = 10.0, delta = 200):
	M = m*UM/1e10
	Rvir = (M/(rhom(1/(1+z))*delta)*3/4/np.pi)**(1/3)
	return 3*GRA*M**2/Rvir/tff(z, delta)/5

def Tvir(m = 1e10, z = 10.0, delta = 200):
	M = m*UM/1e10
	Rvir = (M/(rhom(1/(1+z))*delta)*3/4/np.pi)**(1/3)
	return 3*GRA*M*mmw()*PROTON/Rvir/5/(3*BOL)

import hmf
import hmf.wdm

n_a0 = lambda a: 1.0

ln_M_z = [[20, 3e-5], [19, 8e-5], [18, 2e-4], [17, 7e-4], [16, 2e-3], [15, 5e-3], [14, 1e-2], [13, 3e-2], [12, 6e-2], [11, 0.1], [10, 0.25], [9, 0.35], [8, 0.6], [7, 1.0], [6, 1.5], [5, 2], [4, 2.5], [3, 2.9], [2, 3.0], [1, 2.5], [0, 2.1]]
	
def n_a(a):
	z = 1/a-1
	lref = np.array(ln_M_z).T
	nz = interp1d(lref[0],np.log10(lref[1]))
	return 10**nz(z)

# Radio background
def Jnu_cosmic(z1 = 6, nu = 100, z0 = 20, L = lambda nu: 4e24, n_a = n_a, Om = 0.315, h = 0.6774,mode=0,trec = 100e6):
	if mode==0:
		def integrand(a):
			return L(np.log10(max(min(nu*1e6/a,10*numax()),numin()*10**0.5)))*n_a(a)*SPEEDOFLIGHT*dt_da(a, Om, h)/(UL*1e3/h)**3/4/np.pi
		J = quad(integrand, 1/(1+z0), 1/(1+z1), epsrel = 1e-8)[0]
	elif mode==1:
		def integrand(a):
			return L(np.log10(max(min(nu*1e6/a,10*numax()),numin()*10**0.5)))*SPEEDOFLIGHT * trec * YR * max(derivative(n_a,a,1e-3),0.0) /(UL*1e3/h)**3/4/np.pi #/dt_da(a, Om, h)
		J = quad(integrand, 1/(1+z0), 1/(1+z1), epsrel = 1e-8)[0]
	else:
		J = L(np.log10(max(min(nu*1e6*(1+z1),10*numax()),numin()*10**0.5)))*SPEEDOFLIGHT*n_a(1/(1+z1))/(UL*1e3/h)**3/4/np.pi * trec *YR
	return J*1e23

def Tnu(nu, Jnu):
	return Jnu*SPEEDOFLIGHT**2/2/BOL/(nu*1e6)**2/1e20

def Tnu_sky(nu):
	return 2.725 + 24.1*(nu/310)**-2.599

Ae = Tnu_sky(100)*5*10**3

def Tnu_SKA(nu, dnu = 1.0, t = 1000, sigma = 10):
	return sigma*Tnu_sky(nu)/(2*t*3600*dnu*1e6)**0.5 * 1e3

# Bremsstrahlung
def gff(nu, T):
	lognu = np.log10(nu*HBAR*2*np.pi/T/BOL)
	llognu = np.linspace(-5, 4, 10)
	lgff = [5.8]+[4.5, 3.2, 2.2, 1.7, 1.4, 1.4, 1.5]+[1.45,1.4]
	gff_lognu = interp1d(llognu, lgff, kind='cubic')
	#if lognu<-4:
	#	return -1.3*lognu-0.7
	#elif lognu>2:
	#	a = (1-1/(1+np.exp(3)))
	#	return 1.5  -0.1*(1/(1+np.exp(5-lognu)) -1/(1+np.exp(3)) )/a
	#else:
	#	return gff_lognu(lognu)
	lognu_ = -4*(lognu<-4) + lognu*(-4<=lognu)*(lognu<=2) + 2*(lognu>2)
	a = (1-1/(1+np.exp(3)))
	return (-1.3*lognu-0.7)*(lognu<-4) + gff_lognu(lognu_)*(-4<=lognu)*(lognu<=2) \
			+ (1.5  -0.1*(1/(1+np.exp(5-lognu)) -1/(1+np.exp(3)) )/a)*(lognu>2)

def gff_(nu, T):
	return 11.96*T**0.15*nu**-0.1

def jnu(nu, T, ne, ni, Z = 1):
	enu = 2**5*np.pi*CHARGE**6/3/ELECTRON/SPEEDOFLIGHT**3 * (2*np.pi/3/BOL/ELECTRON)**0.5 / T**0.5 * ne * ni * np.exp(-HBAR*2*np.pi*nu/BOL/T) * gff(nu, T)
	return enu/4/np.pi

def neniT(ne, ni, T):
	return 2**5*np.pi*CHARGE**6/3/ELECTRON/SPEEDOFLIGHT**3 * (2*np.pi/3/BOL/ELECTRON)**0.5 / T**0.5 * ne * ni

def jnu_(nu, T, x):
	return x * np.exp(-HBAR*2*np.pi*nu/BOL/T) * gff(nu, T)/4/np.pi

def Bnu(nu, T):
	x = HBAR*2*np.pi*nu/BOL/T * (HBAR*2*np.pi*nu/BOL/T<100) + 100*(HBAR*2*np.pi*nu/BOL/T>=100)
	return 2*HBAR*2*np.pi*nu**3/SPEEDOFLIGHT**2/(np.exp(x)-1)

# Synchrotron
def jnu_syn(nu, B, ne, p = 2.5, beta = 1.0):
	return beta**2 * 1.3e-17 * ne*(p-1) * B**((p+1)/2) * nu**(-(p-1)/2)

def neB(B, ne, p = 2.5):
	return ne*(p-1) * B**((p+1)/2)

def jnu_syn_(nu, x, p = 2.5, beta = 1.0):
	return beta**2 * 1.3e-17 * x * nu**(-(p-1)/2)

def Snu_syn(nu, B):
	return 2*ELECTRON*nu**2.5 * (2*np.pi*ELECTRON*SPEEDOFLIGHT/B/CHARGE)**0.5

def rat_v(T, v, m = ELECTRON):
	return np.exp(-m*v**2/2/BOL/T)

def Maxwell(T, v, m = ELECTRON):
	A = (m/(2*np.pi*BOL*T))**1.5
	return A*rat_v(T, v, m)

# Mesh
def kernel(r, h):
	q = r/h
	a = 8/np.pi/h**3
	b = ((0<=q)*(q<0.5))*(1+6*q**2*(q-1)) + ((0.5<=q)*(q<1))*2*(1-q)**3
	#if 0<=q<0.5:
	#	b = 1+6*q**2*(q-1)
	#elif 0.5<=q<1:
	#	b = 2*(1-q)**3
	#else:
	#	b = 0
	return a*b

def mesh_3d(x0 = 1750, x1 = 2250, xbin = 100, y0 = 1750, y1 = 2250, ybin = 100, z0 = 1750, z1 = 2250, zbin = 100):
	dx, dy, dz = (x1-x0)/(xbin), (y1-y0)/(ybin), (z1-z0)/(zbin)
	xref = np.linspace(x0+dx/2, x1-dx/2, xbin)
	yref = np.linspace(y0+dy/2, y1-dy/2, ybin)
	zref = np.linspace(z0+dz/2, z1-dz/2, zbin)
	X, Y, Z = np.meshgrid(xref, yref, zref)
	N = xbin*ybin*zbin
	lx = np.reshape(X, N)
	ly = np.reshape(Y, N)
	lz = np.reshape(Z, N)
	lpos = np.array([lx,ly,lz]).T
	d = {}
	d['pos'] = lpos
	d['bins'] = np.array([dx, dy, dz])
	d['bound'] = np.array([[x0,x1],[y0,y1],[z0,z1]])
	d['nb'] = np.array([xbin, ybin, zbin])
	d['mesh'] = [X, Y, Z]
	return d

from coolingf import *
# 0-0 S(0-5), 1-0 Q(1), 1-0 O(3,5), 
# 0-0 S(6-15)
# 1-0 O(2,4,6-8) 
# 1-0 Q(2-7)
# 1-0 S(0-11)
H2_E = np.array([510, 1015, 1682, 2504, 3474, 4586, 6149, 6149, 6956, \
5829, 7197, 8677, 10263, 11940, 13703, 15549, 17458, 19402, 21400, \
5987, 6471, 7584, 8365, 9286, \
6471, 6956, 7584, 8365, 9286, 10341, \
6471, 6956, 7584, 8365, 9286, 10341, 11522, 12817, 14221, 15722, 17311, 18979])

H2_E21 = 1.0/np.array([28.221,17.035,12.279,9.6649,8.0258,6.9091,2.4066,2.8025,3.235, \
6.1088, 5.5115, 5.0529, 4.6947, 4.4096, 4.181, 3.9947, 3.8464, 3.724, 3.625, \
2.6269, 3.0039, 3.5007, 3.8075, 4.1625, \
2.4134, 2.4237, 2.4375, 2.4548, 2.4756, 2.5001, \
2.2235, 2.1218, 2.0338, 1.9576, 1.892, 1.8358, 1.788, 1.748, 1.7147, 1.6877, 1.6665, 1.6504])

H2_g = np.array([5, 21, 9, 33, 13, 45, 9, 9, 21, \
17, 57, 21, 69, 25, 81, 29, 93, 33, 105, \
1, 5, 9, 33, 13, \
5, 21, 9, 33, 13, 45, \
5, 21, 9, 33, 13, 45, 17, 57, 21, 69, 25, 81])

H2_A = np.array([2.94e-11, 4.76e-10,2.76e-9,9.84e-9,2.64e-8,5.88e-8,4.29e-7,4.23e-7,2.09e-7,\
1.14e-7, 2e-7, 3.24e-7, 4.9e-7, 7.03e-7, 9.64e-7, 1.27e-6, 1.62e-6, 2e-6, 2.41e-6, \
8.54e-7, 2.9e-7, 1.5e-7, 1.06e-7, 7.4e-8, \
3.03e-7, 2.78e-7, 2.65e-7, 2.55e-7, 2.45e-7, 2.34e-7,\
2.53e-7, 3.47e-7, 3.98e-7, 4.21e-7, 4.19e-7, 3.96e-7, 3.54e-7, 2.98e-7, 2.34e-7, 1.68e-7, 1.05e-7, 5.3e-8])

def H2_line_dis(T):
	lx2 = np.exp(-H2_E/T)*H2_g
	lf = lx2*H2_E21*H2_A
	return lf/np.sum(lf)

def luminosity_syn(sn, facB = 1.0, facn = 1.0, p_index = 2.5, rep = './', box = [[1900]*3,[2000]*3], Tsh = 1e4, base = 'snapshot', ext = '.hdf5', ncore = 4, X=0.76):
	xh = 4*X/(1+3*X)
	mu0 = 4/(1+3*X)
	ds = yt.load(rep+base+'_'+str(sn).zfill(3)+ext)
	z = ds['Redshift']
	ad = ds.box(box[0],box[1])
	shock = np.array(temp(ad[('PartType0','InternalEnergy')],ad[('PartType0','Primordial HII')]))>Tsh
	nump = ad[('PartType0','Coordinates')][shock].shape[0]
	ln = np.array(ad[('PartType0','Density')][shock].to_equivalent("cm**-3", "number_density",mu=mmw(ad[('PartType0','Primordial HII')][shock])))
	lm = ad[('PartType0','Masses')][shock]
	lrho = ad[('PartType0','Density')][shock]
	lxe = ad[('PartType0','Primordial e-')][shock]
	lB = np.array(magneticf(ad[('PartType0','InternalEnergy')][shock],ad[('PartType0','Density')][shock].to('g/cm**3'),facB))
	print('Maximum field strength: {} [G]'.format(np.max(lB)))
	print('Minimum field strength: {} [G]'.format(np.min(lB)))
	print('Average field strength: {} [G]'.format(np.mean(lB)))
	np_core = int(nump/ncore)
	lpr = [[i*np_core, (i+1)*np_core] for i in range(ncore-1)] + [[(ncore-1)*np_core, nump]]
	print(lpr)
	manager = mp.Manager()
	output = manager.Queue()
	def sess(pr0, pr1):
		out = 0.0
		for i in range(pr0, pr1):
			n = ln[i]*xh
			m = lm[i]
			rho = lrho[i]
			B = lB[i]
			V = np.array((m/rho).to('cm**3'))
			ne = n*lxe[i]*facn
			out += neB(B, ne, p_index)*V
		output.put(out)
	processes = [mp.Process(target=sess, args=(lpr[i][0], lpr[i][1])) for i in range(ncore)]
	for p in processes:
		p.start()
	for p in processes:
		p.join()
	out0 = np.sum([output.get() for p in processes])
	return [z, out0*4*np.pi, p_index]

def luminosity_particle(sn, rep = './', box = [[1900]*3,[2000]*3], nsh = 1e-4, base = 'snapshot', ext = '.hdf5', ncore = 4, X=0.76, nline=42, Tsh = 1e4, nmax = 1.0):
	xh = 4*X/(1+3*X)
	mu0 = 4/(1+3*X)
	ds = yt.load(rep+base+'_'+str(sn).zfill(3)+ext)
	z = ds['Redshift']
	ad = ds.box(box[0],box[1])
	keys = ds.field_list
	tag = np.sum([x[0] == 'PartType3' for x in keys])
	tag_ = np.sum([x[0] == 'PartType4' for x in keys])
	if tag>0:
		Msink = np.array(np.sum(ad[('PartType3','Masses')].to('Msun')))#/ad[('PartType0','Masses')][0])
	else:
		Msink = 0.0
	if tag_>0:
		Msink += np.array(np.sum(ad[('PartType4','Masses')].to('Msun')))
	shock = (ad[('PartType0','Primordial H2')]>nsh) * np.logical_or((np.array(temp(ad[('PartType0','InternalEnergy')],ad[('PartType0','Primordial HII')]))<Tsh),np.array(ad[('PartType0','Density')].to_equivalent("cm**-3", "number_density",mu=mmw(ad[('PartType0','Primordial HII')])))<nmax )
	nump = ad[('PartType0','Coordinates')][shock].shape[0]
	lv = np.array(ad[('PartType0','Velocities')][shock].to('cm/s'))
	ln = np.array(ad[('PartType0','Density')][shock].to_equivalent("cm**-3", "number_density",mu=mmw(ad[('PartType0','Primordial HII')][shock])))
	lm = ad[('PartType0','Masses')][shock]
	lrho = ad[('PartType0','Density')][shock]
	lT = np.array(temp(ad[('PartType0','InternalEnergy')][shock],ad[('PartType0','Primordial HII')][shock]))
	lxH2 = ad[('PartType0','Primordial H2')][shock]
	lxH0 = ad[('PartType0','Primordial HI')][shock]
	lxHD = ad[('PartType0','Primordial HD')][shock]
	np_core = int(nump/ncore)
	lpr = [[i*np_core, (i+1)*np_core] for i in range(ncore-1)] + [[(ncore-1)*np_core, nump]]
	print(lpr)
	manager = mp.Manager()
	output = manager.Queue()
	nu_rest = 1.0#510*BOL/HBAR/2/np.pi
	#print('H2 J=2-0: {} [10^14 Hz]'.format(510*BOL/HBAR/2/np.pi/1e14))
	#print('HD J=1-0: {} [10^14 Hz]'.format(128*BOL/HBAR/2/np.pi/1e14))
	def sess(pr0, pr1):
		ldlh2 = np.zeros(pr1-pr0)
		ldlhd = np.zeros(pr1-pr0)
		lnu_loc = np.zeros(pr1-pr0)
		lnu_scale = np.zeros(pr1-pr0)
		lH2 = np.zeros(nline)
		for i in range(pr0, pr1):
			n = ln[i]*xh
			m = lm[i]
			rho = lrho[i]
			T = lT[i]
			V = np.array((m/rho).to('cm**3'))
			lnu_loc[i-pr0] = nu_rest * lv[i][2]/SPEEDOFLIGHT
			lnu_scale[i-pr0] = nu_rest * (2*BOL*T/PROTON/mu0/SPEEDOFLIGHT**2)**0.5
			#ni = n*ad[('PartType0','Primordial HII')][shock][i]
			#ne = n*ad[('PartType0','Primordial e-')][shock][i]
			nH2 = n*lxH2[i]
			nH0 = n*lxH0[i]
			nHD = n*lxHD[i]*4.3e-5
			ldlh2[i-pr0] = V*LambdaH2(T, nH2, nH0)
			ldlhd[i-pr0] = V*LambdaHD(T, nHD, nH0, n)
			lH2 += ldlh2[i-pr0]*H2_line_dis(T)
			#profile0 = lambda x: profile0(x) + stats.norm.pdf(x, nu_loc, nu_scale)*dlh2
			#profile1 = lambda x: profile1(x) + stats.norm.pdf(x, nu_loc, nu_scale)*dlhd
		output.put([ldlh2, ldlhd, lnu_loc, lnu_scale, lH2])
	processes = [mp.Process(target=sess, args=(lpr[i][0], lpr[i][1])) for i in range(ncore)]
	for p in processes:
		p.start()
	for p in processes:
		p.join()
	out0 = [output.get() for p in processes]
	out_ = np.array([np.hstack([out0[i][j] for i in range(ncore)]) for j in range(4)])
	lH2 = np.sum([out0[i][4] for i in range(ncore)], axis=0)
	out = np.sum(out_[:2],axis=1)
	profile0 = lambda x: np.sum([stats.norm.pdf(x, out_[2][i], out_[3][i])*out_[0][i] for i in range(nump)])/out[0]
	profile1 = lambda x: np.sum([stats.norm.pdf(x, out_[2][i], out_[3][i])*out_[1][i] for i in range(nump)])/out[1]
	print([z, out[0], out[1], Msink])
	print(lH2)
	d = {}
	d['all'] = [z, out[0], out[1], Msink, profile0, profile1]
	d['line'] = lH2
	return d

def grid(mesh, sn, rep = './', box = [[1750]*3,[2250]*3], Tsh = 1e3, mode = 0, base = 'snapshot', ext = '.hdf5', Tmin = 10.0, ncore = 1, X = 0.76, nsh = 1e8):#, h = 0.6774):
	start = time.time()
	ds = yt.load(rep+base+'_'+str(sn).zfill(3)+ext)
	#z = ds['Redshift']
	#posg = mesh['pos']
	dr = mesh['bins']
	ad = ds.box(box[0],box[1])#ds.box([(x, 'kpccm/h') for x in box[0]],[(x, 'kpccm/h') for x in box[1]])#ds.all_data()
	N = mesh['pos'].shape[0]
	shock = (np.array(temp(ad[('PartType0','InternalEnergy')],ad[('PartType0','Primordial HII')]))>Tsh) * (np.array(ad[('PartType0','Density')].to_equivalent("cm**-3", "number_density",mu=mmw(ad[('PartType0','Primordial HII')])))<nsh)
	posp = np.array(ad[('PartType0','Coordinates')].to('kpccm/h'))[shock]
	nump = posp.shape[0]
	print('Number of gas particles (pre): {}'.format(nump))
	ln = np.array(ad[('PartType0','Density')][shock].to_equivalent("cm**-3", "number_density",mu=mmw(ad[('PartType0','Primordial HII')][shock])))
	lm = ad[('PartType0','Masses')][shock]
	lrho = ad[('PartType0','Density')][shock]
	lT = np.array(temp(ad[('PartType0','InternalEnergy')][shock],ad[('PartType0','Primordial HII')][shock]))
	lh = np.array(ad[('PartType0','SmoothingLength')][shock].to('kpccm/h'))
	lxHII = ad[('PartType0','Primordial HII')][shock]
	lxe = ad[('PartType0','Primordial e-')][shock]
	head = [[str(mesh['bound'][x][0]), str(mesh['bound'][x][1]), str(mesh['nb'][x])] for x in range(3)]
	head = head[0]+head[1]+head[2]
	#tcount = 1
	Vc = dr[0]*dr[1]*dr[2]
	low = np.array([mesh['bound'][x][0] for x in range(3)])
	up = np.array([mesh['bound'][x][1] for x in range(3)])
	np_core = int(nump/ncore)
	lpr = [[i*np_core, (i+1)*np_core] for i in range(ncore-1)] + [[(ncore-1)*np_core, nump]]
	print(lpr)
	manager = mp.Manager()
	output = manager.Queue()
	def sess(pr0, pr1):
		mni = np.zeros(mesh['nb'])
		mne = np.zeros(mesh['nb'])
		mT = np.zeros(mesh['nb'])
		meps = np.zeros(mesh['nb'])
		mtag = np.zeros(mesh['nb'])
		mpe = np.zeros(mesh['nb'])
		count = 0
		for i in range(pr0,pr1):
			label = [(posp[i][j]-low[j]-dr[j]/2)/dr[j] for j in range(3)]
			if label[0]<0 or label[1]<0 or label[2]<0:
				continue
			if label[0]>=mesh['nb'][0]-1 or label[1]>=mesh['nb'][1]-1 or label[2]>=mesh['nb'][2]-1:
				continue
			count += 1
			ind = [int(x) for x in label]
			xh = 4*X/(1+3*X+4*X*lxHII[i])
			n = ln[i]*xh
			m = lm[i]
			rho = lrho[i]
			T = lT[i]
			if T<=Tsh:
				print('?')
			h = lh[0]
			V = np.array((m/rho).to('(kpccm/h)**3'))
			ni = n*lxHII[i]
			ne = n*lxe[i]
			for x in [ind[0], ind[0]+1]:
				for y in [ind[1], ind[1]+1]:
					for z in [ind[2], ind[2]+1]:
						posg_ = [x, y, z]
						posg = np.array([(posg_[j]+0.5)*dr[j]+low[j] for j in range(3)])
						dis = np.abs(posp[i] - posg)
						if mode==0:
							side = [dr[j]-dis[j] for j in range(3)]
							weight = side[0]*side[1]*side[2]/Vc
							mni[x][y][z] += ni * weight *V/Vc
							mne[x][y][z] += ne * weight *V/Vc
							mpe[x][y][z] += ne*BOL*T * weight *V/Vc
							mT[x][y][z] += T * weight
							mtag[x][y][z] += weight
							meps[x][y][z] += neniT(ne, ni, T) * weight * V/Vc
							#print(meps[x][y][z])
						else:
							r = np.sum(dis**2)**0.5
							weight = kernel(r,h)*V
							mni[x][y][z] += ni * weight
							mne[x][y][z] += ne * weight
							mpe[x][y][z] += ne*BOL*T * weight
							mT[x][y][z] += T * weight
							mtag[x][y][z] = 1.0
							meps[x][y][z] += neniT(ne, ni, T) * weight
		output.put([count, mni, mne, mT, meps, mpe, mtag])
	processes = [mp.Process(target=sess, args=(lpr[i][0], lpr[i][1])) for i in range(ncore)]
	for p in processes:
		p.start()
	for p in processes:
		p.join()
	mni = np.zeros(mesh['nb'])
	mne = np.zeros(mesh['nb'])
	mT = np.zeros(mesh['nb'])
	meps = np.zeros(mesh['nb'])
	mtag = np.zeros(mesh['nb'])
	mpe = np.zeros(mesh['nb'])
	count = 0
	for p in processes:
		buff = output.get()
		mni += buff[1]
		mne += buff[2]
		mT += buff[3]
		meps += buff[4]
		mpe += buff[5]
		mtag += buff[6]
		count += buff[0]
	lni = np.reshape(mni, N)
	lne = np.reshape(mne, N)
	lT = np.reshape(mT, N)
	leps = np.reshape(meps, N)
	lpe = np.reshape(mpe, N)
	ltag = np.reshape(mtag, N)
	for i in range(N):
		if ltag[i]>0:
			lT[i] = max(lT[i]/ltag[i],Tmin)
		else:
			lT[i] = Tmin
	d = {}
	d['ni'], d['ne'], d['T'] = lni, lne, lT
	d['tag'] = ltag
	d['eps'] = leps
	d['pe'] = lpe
	d['bins'] = mesh['bins']# np.array([dx, dy, dz])
	d['bound'] = mesh['bound']#np.array([[x0,x1],[y0,y1],[z0,z1]])
	d['nb'] = mesh['nb']#np.array([xbin, ybin, zbin])
	d['mesh'] = mesh['mesh']#[X, Y, Z]
	d['ds'] = ds
	nem = len(ltag[ltag==0])
	end = time.time()
	print('Number of gas particles (post): {}'.format(count))
	print('Fraction of empty grids: {}'.format(nem/N))
	print('Time: {} s'.format(end-start))
	totxt(rep+'mesh'+str(mode)+'_'+str(sn)+'_'+str(d['nb'][0])+'.txt',[d['ni'], d['ne'], d['T'], d['eps'], d['pe'], d['tag']],head,1,0)
	return d

def read_grid(sn, nb, rep = './', mode=0, base = 'snapshot', ext = '.hdf5'):
	ds = yt.load(rep+base+'_'+str(sn).zfill(3)+ext)
	l = np.array(retxt(rep+'mesh'+str(mode)+'_'+str(sn)+'_'+str(nb)+'.txt',6,1,0))
	head = restr(rep+'mesh'+str(mode)+'_'+str(sn)+'_'+str(nb)+'.txt')
	x0 = float(head[0][0])
	x1 = float(head[0][1])
	xbin = int(head[0][2])
	y0 = float(head[0][3])
	y1 = float(head[0][4])
	ybin = int(head[0][5])
	z0 = float(head[0][6])
	z1 = float(head[0][7])
	zbin = int(head[0][8])
	dx, dy, dz = (x1-x0)/(xbin), (y1-y0)/(ybin), (z1-z0)/(zbin)
	xref = np.linspace(x0+dx/2, x1-dx/2, xbin)
	yref = np.linspace(y0+dy/2, y1-dy/2, ybin)
	zref = np.linspace(z0+dz/2, z1-dz/2, zbin)
	X, Y, Z = np.meshgrid(xref, yref, zref)
	d = {}
	d['ni'], d['ne'], d['T'] = l[0], l[1], l[2]
	d['eps'] = l[3]
	d['pe'] = l[4]
	d['tag'] = l[5]
	d['bins'] = np.array([dx, dy, dz])
	d['bound'] = np.array([[x0,x1],[y0,y1],[z0,z1]])
	d['nb'] = np.array([xbin, ybin, zbin])
	d['mesh'] = [X, Y, Z]
	d['ds'] = ds
	return d

def numax(T = 1e4):
	return 2.0*BOL*T/HBAR/2/np.pi

def numin(T = 1e4, n = 1e-6, lc = 1000, z = 10, h = 0.6774):
	l = lc*UL/h/(1+z)
	return (4*CHARGE**6/3/ELECTRON/SPEEDOFLIGHT/BOL * (2*np.pi/3/BOL/ELECTRON)**0.5/T**1.5 * n**2 *l)**0.5

def intensity(dg, nu, h = 0.6774):
	z = dg['ds']['Redshift']
	darea = dg['bins'][0]*dg['bins'][1]*(UL/h/(1+z))**2
	ldtau = jnu_(nu, dg['T'], dg['eps'])/Bnu(nu, dg['T']) * dg['bins'][2]*UL/h/(1+z)
	mdtau = np.reshape(ldtau, dg['nb'])
	mpe = np.reshape(dg['pe'], dg['nb'])
	mSnu = Bnu(nu, np.reshape(dg['T'], dg['nb']))
	dx, dy = dg['bins'][0], dg['bins'][1]
	x0, x1 = dg['bound'][0][0], dg['bound'][0][1]
	y0, y1 = dg['bound'][1][0], dg['bound'][1][1]
	xref = np.linspace(x0+dx/2, x1-dx/2, dg['nb'][0])
	yref = np.linspace(y0+dy/2, y1-dy/2, dg['nb'][1])
	X, Y = np.meshgrid(xref, yref)
	Xraymap = np.zeros(X.shape)
	ymap = np.zeros(X.shape)
	for i in range(dg['nb'][0]):
		for j in range(dg['nb'][1]):
			Xraymap[i][j] = np.sum([ np.exp(np.sum(mdtau[i][j][:k])-np.sum(mdtau[i][j])) * mSnu[i][j][k]*mdtau[i][j][k] for k in range(dg['nb'][2])])
			ymap[i][j] = np.sum(mpe[i][j]*SIGMATH/ELECTRON/SPEEDOFLIGHT**2 * dg['bins'][2]*UL/h/(1+z))
	d = {}
	d['X'], d['Y'], d['I'] = X, Y, Xraymap
	d['SZ'] = ymap
	d['L'] = np.sum(Xraymap)*darea*4*np.pi
	return d
			
def SZy(dg, h = 0.6774):
	z = dg['ds']['Redshift']
	mpe = np.reshape(dg['pe'], dg['nb'])
	dx, dy = dg['bins'][0], dg['bins'][1]
	x0, x1 = dg['bound'][0][0], dg['bound'][0][1]
	y0, y1 = dg['bound'][1][0], dg['bound'][1][1]
	xref = np.linspace(x0+dx/2, x1-dx/2, dg['nb'][0])
	yref = np.linspace(y0+dy/2, y1-dy/2, dg['nb'][1])
	X, Y = np.meshgrid(xref, yref)
	ymap = np.zeros(X.shape)
	for i in range(dg['nb'][0]):
		for j in range(dg['nb'][1]):
			ymap[i][j] = np.sum(mpe[i][j]*SIGMATH/ELECTRON/SPEEDOFLIGHT**2 * dg['bins'][2]*UL/h/(1+z))
	d = {}
	d['X'], d['Y'] = X, Y
	d['SZ'] = ymap
	return d

def luminosity(dg, rnu = [np.log10(numax())-2.5, np.log10(numax())+0.5], nb = 100, ncore = 4, h = 0.6774):
	z = dg['ds']['Redshift']
	Vc = (dg['bins'][0] * dg['bins'][1] * dg['bins'][2])*(UL/h/(1+z))**3
	manager = mp.Manager()
	output = manager.Queue()
	lnu0 = np.linspace(10**rnu[0], 10**rnu[1], ncore+1)
	nb_p = int(nb/ncore)
	def sess(i):
		edge = np.linspace(lnu0[i],lnu0[i+1],nb_p+1)
		lnu = (edge[:-1]+edge[1:])/2.0
		dnu = edge[1:]-edge[:-1]
		out = np.sum([Vc*np.sum(jnu_(lnu[x], dg['T'], dg['eps'])*4*np.pi)*dnu[x] for x in range(nb_p)])
		output.put(out)
	processes = [mp.Process(target=sess, args=(i,)) for i in range(ncore)]
	for p in processes:
		p.start()
	for p in processes:
		p.join()
	out = [z, np.sum([output.get() for p in processes])]
	return out
	
def epsnu_bg(nu, z=6, T=2e4, X = 0.76):
	eps = jnu(nu*(1+z)*1e6,T,rhom(1/(1+z))*X/PROTON,rhom(1/(1+z))*X/PROTON)*4*np.pi
	knu = jnu(nu*(1+z)*1e6,T,rhom(1/(1+z))*X/PROTON,rhom(1/(1+z))*X/PROTON)/Bnu(nu*(1+z), T)
	return [eps/(1+z)**3, knu/(1+z)**3]

def Jnu_bg(nu, z0 = 6, z1 = 0, T = 2e4, mode=1):
	if mode==0:
		def dtau(a):
			return epsnu_bg(nu, 1/a-1, T)[1]*dt_da(a)*SPEEDOFLIGHT
		def tau_z(z):
			return quad(dtau, 1/(1+z0), 1/(1+z), epsrel = 1e-8)[0]
		tau0 = tau_z(z1)
		lz = np.linspace(z1,z0,100)
		ltau = np.array([tau_z(x) for x in lz])
		z_tau = interp1d(ltau, lz)
		def integrand(t):
			return Bnu(nu*(1+z_tau(t)), T)*np.exp(-(tau0-t))
		return quad(integrand, 0, tau0, epsrel = 1e-3)[0]*1e23
	else:
		def integrand(a):
			return epsnu_bg(nu, 1/a-1, T)[0]/4/np.pi * dt_da(a)*SPEEDOFLIGHT
		return quad(integrand, 1/(1+z0), 1/(1+z1), epsrel = 1e-8)[0]*1e23

Vbox = (1200*UL/0.6774)**3



