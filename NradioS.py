#from scipy.integrate import quad
#from scipy.interpolate import interp1d
from cosmology import *
import hmf
import sys
#import matplotlib.pyplot as plt
from txt import *

def Nrads_z(lz = np.logspace(-2,np.log10(20),50), nbin = 100, z0 = 30, m1 = 9, m2 = 10, h = 0.6774):
	hmf_ = hmf.MassFunction()
	hmf_.update(n=0.966, sigma_8=0.829,cosmo_params={'Om0':0.315,'H0':67.74},Mmin=m1-1,Mmax=m2+1)
	z0_ = max(z0,np.max(lz))
	z_base = np.linspace(max(0,np.min(lz)-0.1), z0_+0.1, nbin)
	ln_z = np.zeros(nbin)
	for i in range(nbin):
		hmf_.update(z=z_base[i])
		lm = hmf_.m/h
		ln = hmf_.ngtm*h**3
		nm = interp1d(lm,ln)
		ln_z[i] = nm(10**m1)-nm(10**m2)
	logn_z = interp1d(z_base, np.log10(ln_z))
	def integrand(d):
		z = ZD(d)
		n = 10**logn_z(z)
		return 4*np.pi * d**2 * n
	out = np.zeros(lz.shape[0])
	for i in range(lz.shape[0]):
		I = quad(integrand, DZ(lz[i])/1e3/UL, DZ(z0_)/1e3/UL, epsrel = 1e-6)
		out[i] = I[0]
	return [lz, out, [integrand(DZ(x)/1e3/UL) for x in lz]]

if __name__ == "__main__":
	m1 = int(sys.argv[1])
	m2 = int(sys.argv[2])
	z2 = float(sys.argv[3])
	if len(sys.argv)>4:
		z0 = float(sys.argv[4])
	else:
		z0 = 30
	d = Nrads_z(lz = np.logspace(-2,np.log10(z2),50), m1=m1, m2=m2, z0 = z0)
	plt.figure()
	plt.plot(d[0],d[1],label=r'$\int_{d_{C}(z)}^{+\infty}4\pi n(r)r^{2}dr$')
	plt.plot(d[0],d[2],label=r'$4\pi n(d_{C}(z))d_{C}^{2}(z)$',ls='--')
	plt.title(r'$M=10^{'+str(m1)+'}-10^{'+str(m2)+r'}\ M_{\odot}$, $z_{0}='+str(z0)+'$')
	plt.xlabel(r'$z$')
	plt.ylabel(r'Number (per unit comoving distance)')
	plt.yscale('log')
	plt.xlim(0,z2)
	plt.legend()
	plt.tight_layout()
	plt.savefig('Nrads_z_'+str(m1)+'_'+str(m2)+'.pdf')
	plt.show()
