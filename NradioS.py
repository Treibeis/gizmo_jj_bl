#from scipy.integrate import quad
#from scipy.interpolate import interp1d
from cosmology import *
import hmf
import sys
import matplotlib.pyplot as plt
#from txt import *
import argparse

def Nrads_z(lz = np.logspace(-2,np.log10(20),50), nbin = 100, z0 = 30, m1 = 9, m2 = 10, trec = 100, h = 0.6774):
	hmf_ = hmf.MassFunction()
	hmf_.update(n=0.966, sigma_8=0.829,cosmo_params={'Om0':0.315,'H0':67.74},Mmin=m1-1,Mmax=m2+1)
	z0_ = max(z0,np.max(lz))
	#z_base = np.linspace(max(0,np.min(lz)-0.1), z0_+0.1, nbin)
	z_base = [ZT(x) for x in np.log10(np.linspace(TZ(max(0,np.min(lz)-0.1)), TZ(z0_+0.1), nbin)/1e9/YR)]
	ln_z = np.zeros(nbin)
	for i in range(nbin):
		hmf_.update(z=z_base[i])
		lm = hmf_.m/h
		ln = hmf_.ngtm*h**3
		nm = interp1d(lm,ln)
		nf = nm(10**m1)-nm(10**m2)
		zi = ZT(np.log10( max((TZ(z_base[i])-trec*1e6*YR)/1e9/YR, 0.2)))
		hmf_.update(z=zi)
		lm = hmf_.m/h
		ln = hmf_.ngtm*h**3
		nm = interp1d(lm,ln)
		ni = nm(10**m1)-nm(10**m2)
		if zi>z_base[i]:
			ln_z[i] = max(nf-ni,0)
		else:
			ln_z[i] = 0.0
	n_z = interp1d(z_base, ln_z)
	def integrand(d):
		z = ZD(d)
		n = n_z(z)
		return 4*np.pi * d**2 * n
	out = np.zeros(lz.shape[0])
	for i in range(lz.shape[0]):
		I = quad(integrand, DZ(lz[i])/1e3/UL, DZ(z0_)/1e3/UL, epsrel = 1e-6)
		out[i] = I[0]
	return [lz, out, [integrand(DZ(x)/1e3/UL) for x in lz]]

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Number of DM halos in the sky')
	parser.add_argument('--m1', help='log10(m1 [Msun])', type=float, required=False)
	parser.add_argument('--m2', help='log10(m2 [Msun])', type=float, required=False)
	parser.add_argument('--z1', help='z1', type=float, required=False)
	parser.add_argument('--z2', help='z2', type=float, required=False)
	parser.add_argument('--z0', help='start redshift', type=float, required=False)
	parser.add_argument('--t', help='time scale [Myr]', type=float, required=False)
	args = parser.parse_args()
	args = vars(args)
	if args['m1']==None:
		m1 = 9.0
	else:
		m1 = args['m1']
	if args['m2']==None:
		m2 = 10.0
	else:
		m2 = args['m2']
	if args['z1']==None:
		z1 = 0.01
	else:
		z1 = args['z1']
	if args['z2']==None:
		z2 = 20.0
	else:
		z2 = args['z2']
	if args['z0']==None:
		z0 = 30.0
	else:
		z0 = args['z0']
	if args['t']==None:
		trec = 100.0
	else:
		trec = args['t']
	d = Nrads_z(lz = np.logspace(np.log10(z1),np.log10(z2),50), m1=m1, m2=m2, z0 = z0, trec=trec)
	plt.figure()
	plt.plot(d[0],d[1],label=r'$\int_{d_{C}(z)}^{+\infty}4\pi n_{\tau}(r)r^{2}dr$')
	plt.plot(d[0],d[2],label=r'$4\pi n_{\tau}(d_{C}(z))d_{C}^{2}(z)$',ls='--')
	plt.title(r'$M=10^{'+str(m1)+'}-10^{'+str(m2)+r'}\ M_{\odot}$, $z_{0}='+str(z0)+r'$, $\tau='+str(trec)+'\ \mathrm{Myr}$')
	plt.xlabel(r'$z$')
	plt.ylabel(r'Number (per unit comoving distance)')
	plt.yscale('log')
	plt.xlim(z1,z2)
	plt.legend()
	plt.tight_layout()
	plt.savefig('Nrads_z_'+str(m1)+'_'+str(m2)+'.pdf')
	plt.show()
