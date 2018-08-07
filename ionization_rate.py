from txt import *
from scipy.interpolate import interp1d

if __name__ == "__main__":
	zmin, zmax = 1, 10
	zbins = 100
	zbase = np.linspace(zmin, zmax, zbins)

	tab = retxt('TREECOOL',6,0,0)
	llogz = tab[0]
	llogzeta = np.log10(tab[1])
	logzeta = interp1d(llogz, llogzeta)
	zeta = lambda z: 10**logzeta(np.log10(z))
	lzeta = zeta(zbase)
	print('zeta(z=10.23): {}, zeta(z=9.84): {}'.format(zeta(10.23),zeta(9.84)))

	plt.figure()
	plt.plot(zbase,lzeta)
	#plt.yscale('log')
	plt.xlabel(r'$z$')
	plt.ylabel(r'$\zeta\ [\mathrm{s^{-1}}]$')
	plt.tight_layout()
	plt.savefig('zeta_z.pdf')
	plt.show()
