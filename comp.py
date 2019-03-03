from radio import *

lh = [1, 3, 7]
lls = ['-', '--', '-.']

def main(lh, ref_ind = 16, tcore = 10, bins = 100, low=1750, up=2250):
	lrep = ['halo'+str(x)+'_jj/' for x in lh]
	lout0 = [np.array(retxt(rep+'luminosity_z_'+str(bins)+'_'+lmodel[1]+'.txt',2,1,0)) for rep in lrep]
	lout1 = [np.array(retxt(rep+'luminosity_z_'+str(bins)+'_'+lmodel[0]+'.txt',2,1,0)) for rep in lrep]
	llu0 = [np.array(retxt(rep+'luminosityH2_z_'+str(low)+'_'+str(up)+'_'+lmodel[1]+'.txt',4,1,0)) for rep in lrep]
	llu1 = [np.array(retxt(rep+'luminosityH2_z_'+str(low)+'_'+str(up)+'_'+lmodel[0]+'.txt',4,1,0)) for rep in lrep]
	lltot0 = [np.array(retxt(rep+'Ltot_z_'+str(bins)+'_'+lmodel[1]+'.txt',5,1,0)) for rep in lrep]
	lltot1 = [np.array(retxt(rep+'Ltot_z_'+str(bins)+'_'+lmodel[0]+'.txt',5,1,0)) for rep in lrep]
	llz0 = []
	llz1 = []
	llff0 = []
	llff1 = []
	Mref0 = lltot0[0][2][ref_ind]
	Mref1 = lltot1[0][2][ref_ind]
	for i in range(len(lh)):
		out0, out1 = lout0[i], lout1[i]
		lu0, lu1 = llu0[i], llu1[i]
		ltot0, ltot1 = lltot0[i], lltot1[i]
		lt0 = np.array([TZ(x)/YR/1e6 for x in lu0[0]])
		lt1 = np.array([TZ(x)/YR/1e6 for x in lu1[0]])
		Ms_t0 = interp1d(lt0, lu0[3])
		Ms_t1 = interp1d(lt1, lu1[3])
		lt0[1:] = lt0[1:]-tcore
		lt1[1:] = lt1[1:]-tcore
		luc0 = lu0[3]-Ms_t0(lt0)
		luc1 = lu1[3]-Ms_t1(lt1)
		lHII0 = epsilon_ff*luc0*UM/(1e10*mmw()*100*PROTON)
		lHII1 = epsilon_ff*luc1*UM/(1e10*mmw()*100*PROTON)
		llz0.append(lu0[0])
		llz1.append(lu1[0])
		llff0.append((out0[1] + lHII0)*ltot0[2][ref_ind]/Mref0)
		llff1.append((out1[1] + lHII1)*ltot1[2][ref_ind]/Mref1)
	plt.figure(figsize=(12,5))
	plt.subplot(121)
	a = [plt.plot(llz0[i][llff0[i]>0], llff0[i][llff0[i]>0], ls = lls[i], label='Halo_'+str(lh[i])) for i in range(len(lh))]
	plt.xlabel(r'$z$')
	plt.ylabel(r'$\hat{L}_{\mathrm{ff}}\ [\mathrm{erg\ s^{-1}}]$')
	plt.yscale('log')
	plt.legend()
	plt.title('WDM_3_KeV')
	plt.subplot(122)
	a = [plt.plot(llz1[i][llff1[i]>0], llff1[i][llff1[i]>0], ls = lls[i], label='Halo_'+str(lh[i])) for i in range(len(lh))]
	plt.xlabel(r'$z$')
	plt.ylabel(r'$\hat{L}_{\mathrm{ff}}\ [\mathrm{erg\ s^{-1}}]$')
	plt.yscale('log')
	plt.legend()
	plt.title('CDM')
	plt.tight_layout()
	plt.savefig('Lff_var.pdf')
	plt.close()
	for i in range(len(llz0[0])):
		if np.sum([int(x[i]>0) for x in llff0])>=3:
			ind0 = i+1
			break
	for i in range(len(llz1[0])):
		if np.sum([int(x[i]>0) for x in llff1])>=3:
			ind1 = i+1
			break
	lrat0 = np.array([x[ind0:ref_ind+1]/llff0[0][ind0:ref_ind+1] for x in llff0]).T
	lrat1 = np.array([x[ind1:ref_ind+1]/llff1[0][ind1:ref_ind+1] for x in llff1]).T
	lflu0 = [np.std(x) for x in lrat0]
	lflu1 = [np.std(x) for x in lrat1]
	print(lflu0, np.min(lflu0), np.max(lflu0))
	print(lflu1, np.min(lflu1), np.max(lflu1))

main(lh)

lpos = [[0.71811605, 0.39721545, 0.2059604785], [0.5543034, 0.45781925, 0.851971], [0.5858542, 0.5209606, 0.89887995]]
for i in range(3):
	for j in range(i+1,3):
		p0 = lpos[i]
		p1 = lpos[j]
		x = np.abs(p0[0]-p1[0])
		y = np.abs(p0[1]-p1[1])
		z = np.abs(p0[2]-p1[2])
		if x>0.5:
			x = 1.-x
		if y>0.5:
			y = 1.-y
		if z>0.5:
			z = 1.-z
		d = 4.*(x**2+y**2+z**2)**0.5
		print('d({}, {}) = {:.3f} h^-1 Mpc'.format(i,j,d))


