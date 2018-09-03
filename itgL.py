from radio import *

if __name__ == "__main__":
	tag = 1
	sca = 1
	sfdbk = 1

	ncore = 6
	nline = 42
	#rep0 = 'halo1_jj/'
	rep0 = 'halo1/'
	ldir = ['NL4_zoom_wdm/'+rep0, 'NL4_zoom_cdm/'+rep0]
	#ldir = ['halo1_wdm/','halo1_cdm/']
	Tsh = 1e4
	mode = int(sys.argv[1])
	bins = int(sys.argv[2])
	if len(sys.argv)>4:
		low, up = int(sys.argv[3]), int(sys.argv[4])
	else:
		low, up = 1900, 2000

	sn0 = 23#5
	sn1 = 23#5

	if tag==0:
		out0 = []
		out1 = []
		dlu0, dlu1 = [], []
		dltot0, dltot1 = [], []
		for sn in range(0,26):
			if mode==0:
				mesh = mesh_3d(low, up, bins, low, up, bins, low, up, bins)
				if sn<=sn0:
					d0 = grid(mesh,sn,box =[[low]*3,[up]*3],ncore=ncore,Tsh=Tsh,rep=ldir[0])
				if sn<=sn1:
					d1 = grid(mesh,sn,box =[[low]*3,[up]*3],ncore=ncore,Tsh=Tsh,rep=ldir[1])
			else:
				if sn<=sn0:
					d0 = read_grid(sn,nb = bins,rep=ldir[0])
				if sn<=sn1:
					d1 = read_grid(sn,nb = bins,rep=ldir[1])
			if sn<=sn0:
				out0.append(luminosity(d0,ncore=ncore))
				lu0_temp = luminosity_particle(sn,ldir[0],[[low]*3,[up]*3],ncore=ncore)
				dlu0.append(lu0_temp)
			if sn<=sn1:
				out1.append(luminosity(d1,ncore=ncore))
				lu1_temp = luminosity_particle(sn,ldir[1],[[low]*3,[up]*3],ncore=ncore)		
				dlu1.append(lu1_temp)
		out0 = np.array(out0).T
		out1 = np.array(out1).T
		lu0 = np.array([x['all'][:4] for x in dlu0]).T
		lu1 = np.array([x['all'][:4] for x in dlu1]).T
		totxt(rep0+'luminosity_z_'+str(bins)+'_'+lmodel[1]+'.txt',out0,['z', 'L(WDM)'],1,0)
		totxt(rep0+'luminosity_z_'+str(bins)+'_'+lmodel[0]+'.txt',out1,['z', 'L(CDM)'],1,0)
		totxt(rep0+'luminosityH2_z_'+str(low)+'_'+str(up)+'_'+lmodel[1]+'.txt',lu0,['z', 'L(WDM)_H2', 'L(WDM)_HD', 'ns(WDM)'],1,0)
		totxt(rep0+'luminosityH2_z_'+str(low)+'_'+str(up)+'_'+lmodel[0]+'.txt',lu1,['z', 'L(CDM)_H2', 'L(CDM)_HD', 'ns(CDM)'],1,0)
		lH20 = np.array([x['line'] for x in dlu0]).T
		lH21 = np.array([x['line'] for x in dlu1]).T
		totxt(rep0+'luminosityH2_line_z_'+str(low)+'_'+str(up)+'_'+lmodel[1]+'.txt',[lu0[0]]+[x for x in lH20],['z']+['0-0 S('+str(i)+')' for i in range(6)]+['1-0 Q(1)', '1-0 O(3)', '1-0 O(5)'],1,0)
		totxt(rep0+'luminosityH2_line_z_'+str(low)+'_'+str(up)+'_'+lmodel[0]+'.txt',[lu1[0]]+[x for x in lH21],['z']+['0-0 S('+str(i)+')' for i in range(6)]+['1-0 Q(1)', '1-0 O(3)', '1-0 O(5)'],1,0)
	else:
		out0 = np.array(retxt(rep0+'luminosity_z_'+str(bins)+'_'+lmodel[1]+'.txt',2,1,0))
		out1 = np.array(retxt(rep0+'luminosity_z_'+str(bins)+'_'+lmodel[0]+'.txt',2,1,0))
	lu0 = np.array(retxt(rep0+'luminosityH2_z_'+str(low)+'_'+str(up)+'_'+lmodel[1]+'.txt',4,1,0))
	lu1 = np.array(retxt(rep0+'luminosityH2_z_'+str(low)+'_'+str(up)+'_'+lmodel[0]+'.txt',4,1,0))
	lH20 = np.array(retxt(rep0+'luminosityH2_line_z_'+str(low)+'_'+str(up)+'_'+lmodel[1]+'.txt',nline+1,1,0))
	lH21 = np.array(retxt(rep0+'luminosityH2_line_z_'+str(low)+'_'+str(up)+'_'+lmodel[0]+'.txt',nline+1,1,0))
	#print(lH20)

	plt.figure()
	plt.plot(lu0[0][lu0[1]>0],lH20[2][lu0[1]>0],label='0-0 S(1), '+lmodel_[1],marker='o')
	plt.plot(lu1[0][lu1[1]>0],lH21[2][lu1[1]>0],label='0-0 S(1), '+lmodel_[0],ls='--',marker='o')
	plt.plot(lu0[0][lu0[1]>0],lH20[4][lu0[1]>0],label='0-0 S(3), '+lmodel_[1],marker='^')
	plt.plot(lu1[0][lu1[1]>0],lH21[4][lu1[1]>0],label='0-0 S(3), '+lmodel_[0],ls='--',marker='^')
	plt.xlabel(r'$z$')
	plt.xlim(min(lu0[0][lu0[1]>0])-0.1,max(lu1[0][lu1[1]>0])+0.1)
	plt.ylabel(r'$L_{\mathrm{H_{2}}}\ [\mathrm{erg\ s^{-1}}]$')
	if sca!=0:
		plt.yscale('log')
		#plt.xscale('log')
	plt.legend()
	plt.tight_layout()
	if sca==0:
		plt.savefig(rep0+'LH2_line_z_'+str(bins)+'.pdf')
	else:
		plt.savefig(rep0+'logLH2_line_z_'+str(bins)+'.pdf')
	#plt.show()

	lz0 = lu0[0][lu0[1]>0]
	lz1 = lu1[0][lu1[1]>0]
	lfs10 = [lH20[2][lu0[1]>0][i]/(DZ(lz0[i])*(1+lz0[i]))**2/4/np.pi/1e3 for i in range(len(lz0))]
	lfs11 = [lH21[2][lu1[1]>0][i]/(DZ(lz1[i])*(1+lz1[i]))**2/4/np.pi/1e3 for i in range(len(lz1))]
	lfs30 = [lH20[4][lu0[1]>0][i]/(DZ(lz0[i])*(1+lz0[i]))**2/4/np.pi/1e3 for i in range(len(lz0))]
	lfs31 = [lH21[4][lu1[1]>0][i]/(DZ(lz1[i])*(1+lz1[i]))**2/4/np.pi/1e3 for i in range(len(lz1))]
	lfs50 = [lH20[6][lu0[1]>0][i]/(DZ(lz0[i])*(1+lz0[i]))**2/4/np.pi/1e3 for i in range(len(lz0))]
	lfs51 = [lH21[6][lu1[1]>0][i]/(DZ(lz1[i])*(1+lz1[i]))**2/4/np.pi/1e3 for i in range(len(lz1))]
	#lfs70 = [lH20[17][lu0[1]>0][i]/(DZ(lz0[i])*(1+lz0[i]))**2/4/np.pi/1e3 for i in range(len(lz0))]
	#lfs71 = [lH21[17][lu1[1]>0][i]/(DZ(lz1[i])*(1+lz1[i]))**2/4/np.pi/1e3 for i in range(len(lz1))]
	plt.figure()
	plt.plot(lz0,lfs10,label='0-0 S(1), '+lmodel_[1],marker='o')
	plt.plot(lz0,lfs30,label='0-0 S(3), '+lmodel_[1],marker='^')
	plt.plot(lz0,lfs50,label='0-0 S(5), '+lmodel_[1],marker='*')
	plt.plot(lz1,lfs11,label='0-0 S(1), '+lmodel_[0],ls='--',marker='o')
	plt.plot(lz1,lfs31,label='0-0 S(3), '+lmodel_[0],ls='--',marker='^')
	plt.plot(lz1,lfs51,label='0-0 S(5), '+lmodel_[0],ls='--',marker='*')
	#plt.plot(lz0,lfs70,label='0-0 S(13), '+lmodel[1],marker='+')
	#plt.plot(lz1,lfs71,label='0-0 S(13), '+lmodel[0],ls='--',marker='+')
	plt.xlabel(r'$z$')
	plt.xlim(min(lu0[0][lu0[1]>0])-0.1,max(lu1[0][lu1[1]>0])+0.1)
	plt.ylabel(r'$F_{\mathrm{H_{2}}}\ [\mathrm{W\ m^{-2}}]$')
	#plt.title(r'$M_{\mathrm{vir}}=6.9\ (7.67)\times 10^{9}\ M_{\odot}$ in WDM (CDM) cosmology,'+'\n with stellar feedbacks',size=12)
	if sca!=0:
		plt.yscale('log')
		#plt.xscale('log')
	plt.legend()
	plt.tight_layout()
	if sca==0:
		plt.savefig(rep0+'FH2_line_z_'+str(bins)+'.pdf')
	else:
		plt.savefig(rep0+'logFH2_line_z_'+str(bins)+'.pdf')

	#lflux0 = [llu0[i]/(DZ(lz0[i])*(1+lz0[i]))**2/4/np.pi/1e3 for i in range(len(lz0))]
	#lflux1 = [llu1[i]/(DZ(lz1[i])*(1+lz1[i]))**2/4/np.pi/1e3 for i in range(len(lz1))]

	plt.figure()
	plt.plot(lu0[0][lu0[1]>0],lH20[4][lu0[1]>0],label='0-0 S(3), '+lmodel[1],marker='^')
	#plt.plot(lu1[0][lu1[1]>0],lH21[4][lu1[1]>0],label='0-0 S(3), '+lmodel[0],ls='--',marker='^')
	plt.plot(lu0[0][lu0[1]>0],lH20[6][lu0[1]>0],label='0-0 S(5), '+lmodel[1],marker='.')
	#plt.plot(lu1[0][lu1[1]>0],lH21[6][lu1[1]>0],label='0-0 S(5), '+lmodel[0],ls='--',marker='+')
	plt.plot(lu0[0][lu0[1]>0],lH20[7][lu0[1]>0],label='1-0 Q(1), '+lmodel[1],marker='x')
	#plt.plot(lu1[0][lu1[1]>0],lH21[7][lu1[1]>0],label='1-0 Q(1), '+lmodel[0],ls='--',marker='x')
	plt.plot(lu0[0][lu0[1]>0],lH20[8][lu0[1]>0],label='1-0 O(3), '+lmodel[1],marker='*')
	#plt.plot(lu1[0][lu1[1]>0],lH21[8][lu1[1]>0],label='1-0 O(3), '+lmodel[0],ls='--',marker='*')
	plt.plot(lu0[0][lu0[1]>0],lH20[9][lu0[1]>0],label='1-0 O(5), '+lmodel[1],marker='o')
	plt.xlabel(r'$z$')
	plt.xlim(min(lu0[0][lu0[1]>0])-0.1,max(lu1[0][lu1[1]>0])+0.1)
	plt.ylabel(r'$L_{\mathrm{H_{2}}}\ [\mathrm{erg\ s^{-1}}]$')
	if sca!=0:
		plt.yscale('log')
		#plt.xscale('log')
	plt.legend()
	plt.tight_layout()
	if sca==0:
		plt.savefig(rep0+'LH2_line_WDM_z_'+str(bins)+'.pdf')
	else:
		plt.savefig(rep0+'logLH2_line_WDM_z_'+str(bins)+'.pdf')

	plt.figure()
	#plt.plot(lu0[0][lu0[1]>0],lH20[4][lu0[1]>0],label='0-0 S(3), '+lmodel[1],marker='^')
	plt.plot(lu1[0][lu1[1]>0],lH21[4][lu1[1]>0],label='0-0 S(3), '+lmodel[0],ls='--',marker='^')
	#plt.plot(lu0[0][lu0[1]>0],lH20[6][lu0[1]>0],label='0-0 S(5), '+lmodel[1],marker='+')
	plt.plot(lu1[0][lu1[1]>0],lH21[6][lu1[1]>0],label='0-0 S(5), '+lmodel[0],ls='--',marker='.')
	#plt.plot(lu0[0][lu0[1]>0],lH20[7][lu0[1]>0],label='1-0 Q(1), '+lmodel[1],marker='x')
	plt.plot(lu1[0][lu1[1]>0],lH21[7][lu1[1]>0],label='1-0 Q(1), '+lmodel[0],ls='--',marker='x')
	#plt.plot(lu0[0][lu0[1]>0],lH20[8][lu0[1]>0],label='1-0 O(3), '+lmodel[1],marker='*')
	plt.plot(lu1[0][lu1[1]>0],lH21[8][lu1[1]>0],label='1-0 O(3), '+lmodel[0],ls='--',marker='*')
	plt.plot(lu1[0][lu1[1]>0],lH21[9][lu1[1]>0],label='1-0 O(5), '+lmodel[0],ls='--',marker='o')
	plt.xlabel(r'$z$')
	plt.xlim(min(lu0[0][lu0[1]>0])-0.1,max(lu1[0][lu1[1]>0])+0.1)
	plt.ylabel(r'$L_{\mathrm{H_{2}}}\ [\mathrm{erg\ s^{-1}}]$')
	if sca!=0:
		plt.yscale('log')
		#plt.xscale('log')
	plt.legend()
	plt.tight_layout()
	if sca==0:
		plt.savefig(rep0+'LH2_line_CDM_z_'+str(bins)+'.pdf')
	else:
		plt.savefig(rep0+'logLH2_line_CDM_z_'+str(bins)+'.pdf')

	lt0 = np.array([TZ(x)/YR/1e6 for x in lu0[0]])
	lt1 = np.array([TZ(x)/YR/1e6 for x in lu1[0]])
	Ms_t0 = interp1d(lt0, lu0[3])
	Ms_t1 = interp1d(lt1, lu1[3])
	tcore = 0.1
	lt0[1:] = lt0[1:]-tcore
	lt1[1:] = lt1[1:]-tcore
	luc0 = lu0[3]-Ms_t0(lt0)
	luc1 = lu1[3]-Ms_t1(lt1)

	lz0 = lu0[0][lu0[1]>0]
	lz1 = lu1[0][lu1[1]>0]
	llu0 = (lu0[1]+luc0*0.1*5e33/10)[lu0[1]>0]
	llu1 = (lu1[1]+luc1*0.1*5e33/10)[lu1[1]>0]
	lflux0 = [llu0[i]/(DZ(lz0[i])*(1+lz0[i]))**2/4/np.pi/1e3 for i in range(len(lz0))]
	lflux1 = [llu1[i]/(DZ(lz1[i])*(1+lz1[i]))**2/4/np.pi/1e3 for i in range(len(lz1))]
	plt.figure()
	plt.plot(lz0,lflux0,label='diffuse+core, '+lmodel_[1],marker='*')
	plt.plot(lz1,lflux1,label='diffuse+core, '+lmodel_[0],ls='--',marker='*')
	#plt.plot(lu0[0][lu0[3]>0],lu0[3][lu0[3]>0]*0.05*5e33/10,label='core, '+lmodel[1])
	#plt.plot(lu1[0][lu1[3]>0],lu1[3][lu1[3]>0]*0.05*5e33/10,label='core ($\epsilon=0.05$, $M_{*}=10\ M_{\odot}$), '+lmodel[0],ls='--')
	plt.xlabel(r'$z$')
	plt.xlim(min(lz0)-0.1,max(lz1)+0.1)
	plt.ylabel(r'$F_{\mathrm{H_{2}}}\ [\mathrm{W\ m^{-2}}]$')
	if sca!=0:
		plt.yscale('log')
		#plt.xscale('log')
	plt.legend()
	plt.tight_layout()
	if sca==0:
		plt.savefig(rep0+'FH2_z_'+str(bins)+'.pdf')
	else:
		plt.savefig(rep0+'logFH2_z_'+str(bins)+'.pdf')
	print('H2 flux: {} (CDM), {} (WDM) [W m^-2]'.format(max(lflux1),max(lflux0)))
	print('H2 flux at z = {}: {} (CDM), {} (WDM) [W m^-2]'.format(lu0[0][18],lflux1[-8],lflux0[-8]))

	if tag!=0:
		ref_ind = 20
		Mv_ = 7e9
		ltot0 = np.array(retxt(rep0+'Ltot_z_'+str(bins)+'_'+lmodel[1]+'.txt',4,1,0))
		ltot1 = np.array(retxt(rep0+'Ltot_z_'+str(bins)+'_'+lmodel[0]+'.txt',4,1,0))
		lMvir0 = np.array([Mv_*ltot0[2][i]/ltot0[2][ref_ind] for i in range(ltot0.shape[1])])
		lMvir1 = np.array([Mv_*ltot1[2][i]/ltot1[2][ref_ind] for i in range(ltot1.shape[1])])
		lvir0 = np.array([Lvir(lMvir0[i],ltot0[0][i]) for i in range(ltot0.shape[1])])
		lvir1 = np.array([Lvir(lMvir1[i],ltot1[0][i]) for i in range(ltot1.shape[1])])

	plt.figure()
	plt.plot(out0[0][out0[1]>0],out0[1][out0[1]>0],label=r'$L_{\mathrm{ff}}$, '+lmodel_[1],marker='^')
	plt.plot(out1[0][out1[1]>0],out1[1][out1[1]>0],label=r'$L_{\mathrm{ff}}$, '+lmodel_[0],ls='--',marker='^')
	plt.plot(lu0[0][lu0[1]>0],lu0[1][lu0[1]>0],label=r'$L_{\mathrm{H_{2}}}^{\mathrm{D}}$, '+lmodel_[1],marker='o')
	plt.plot(lu1[0][lu1[1]>0],lu1[1][lu1[1]>0],label=r'$L_{\mathrm{H_{2}}}^{\mathrm{D}}$, '+lmodel_[0],ls='--',marker='o')
	plt.plot(lu0[0][lu0[3]>0],luc0[lu0[3]>0]*0.1*5e33/10,label=r'$L_{\mathrm{H_{2}}}^{\mathrm{C}}$, '+lmodel_[1],marker='.')
	plt.plot(lu1[0][lu1[3]>0],luc1[lu1[3]>0]*0.1*5e33/10,label=r'$L_{\mathrm{H_{2}}}^{\mathrm{C}}$, '+lmodel_[0],ls='--',marker='.')# ($\epsilon=0.05$, $M_{*}=10\ M_{\odot}$)
	if tag!=0:
		plt.plot(ltot0[0][ltot0[1]>0],ltot0[1][ltot0[1]>0],label=r'$L_{\mathrm{tot}}$, '+lmodel_[1],marker='*')#,lw=1)
		plt.plot(ltot1[0][ltot1[1]>0],ltot1[1][ltot1[1]>0],label=r'$L_{\mathrm{tot}}$, '+lmodel_[0],marker='*',ls='--')#,lw=1)
		plt.plot(ltot0[0][ltot0[2]>0],lvir0[ltot0[2]>0],label=r'$L_{\mathrm{vir}}$, '+lmodel_[1],marker='x',lw=1)
		plt.plot(ltot1[0][ltot1[2]>0],lvir1[ltot1[2]>0],label=r'$L_{\mathrm{vir}}$, '+lmodel_[0],marker='x',ls='--',lw=1)
		plt.plot(ltot0[0][ltot0[2]>0],3.2e4*lMvir0[ltot0[2]>0]*2e33,label=r'$L_{\mathrm{Edd}}$, '+lmodel_[1],lw=3)
		plt.plot(ltot1[0][ltot1[2]>0],3.2e4*lMvir1[ltot1[2]>0]*2e33,label=r'$L_{\mathrm{Edd}}$, '+lmodel_[0],ls='--',lw=3)
	plt.xlabel(r'$z$')
	#plt.xlim(min(lu0[0][lu0[1]>0])-0.1,max(lu1[0][lu1[1]>0])+0.1)
	plt.xlim(min(lu0[0][lu0[1]>0])-0.1,30)
	plt.ylabel(r'$L\ [\mathrm{erg\ s^{-1}}]$')
	if sca!=0:
		plt.yscale('log')
		#plt.xscale('log')
	plt.legend()
	plt.tight_layout()
	if sca==0:
		plt.savefig(rep0+'Ltot_z_'+str(bins)+'.pdf')
	else:
		plt.savefig(rep0+'logLtot_z_'+str(bins)+'.pdf')

	plt.figure()
	plt.plot(lu0[0][lu0[2]>0],lu0[2][lu0[2]>0],label=lmodel_[1],marker='^')
	plt.plot(lu1[0][lu1[2]>0],lu1[2][lu1[2]>0],label=lmodel_[0],ls='--',marker='^')
	plt.xlabel(r'$z$')
	plt.xlim(min(lu0[0][lu0[2]>0])-0.1,max(lu1[0][lu1[2]>0])+0.1)
	plt.ylabel(r'$L_{\mathrm{HD}}\ [\mathrm{erg\ s^{-1}}]$')
	if sca!=0:
		plt.yscale('log')
		#plt.xscale('log')
	plt.legend()
	plt.tight_layout()
	if sca==0:
		plt.savefig(rep0+'LHD_z_'+str(bins)+'.pdf')
	else:
		plt.savefig(rep0+'logLHD_z_'+str(bins)+'.pdf')

	plt.figure()
	plt.plot(lu0[0][lu0[3]>0],lu0[3][lu0[3]>0],label=lmodel_[1],marker='o')
	plt.plot(lu1[0][lu1[3]>0],lu1[3][lu1[3]>0],label=lmodel_[0],ls='--',marker='o')
	plt.xlabel(r'$z$')
	plt.ylabel(r'$M_{\mathrm{sink}}\ [M_{\odot}]$')
	if sca!=0:
		plt.yscale('log')
		#plt.xscale('log')
	plt.legend()
	plt.tight_layout()
	if sca==0:
		plt.savefig(rep0+'Msink_z_'+str(bins)+'.pdf')
	else:
		plt.savefig(rep0+'logMsink_z_'+str(bins)+'.pdf')
	#plt.show()
	
	tcore = 10
	lt0[1:] = lt0[1:]-tcore
	lt1[1:] = lt1[1:]-tcore
	luc0 = lu0[3]-Ms_t0(lt0)
	luc1 = lu1[3]-Ms_t1(lt1)

	lz0 = lu0[0][lu0[3]>0]
	lz1 = lu1[0][lu1[3]>0]
	lHII0 = epsilon_ff*luc0*UM/(1e10*mmw()*100*PROTON)
	lHII1 = epsilon_ff*luc1*UM/(1e10*mmw()*100*PROTON)
	print('epsilon_ff: {} [erg s^-1]'.format(epsilon_ff))

	plt.figure()
	plt.plot(out0[0][out0[1]>0],out0[1][out0[1]>0],label=lmodel_[1],marker='^')
	plt.plot(out1[0][out1[1]>0],out1[1][out1[1]>0],label=lmodel_[0],ls='--',marker='^')
	if sfdbk!=0:
		plt.plot(lz0, lHII0[lu0[3]>0],label=r'Unresolved, '+lmodel_[1],marker='o')
		plt.plot(lz1, lHII1[lu1[3]>0],label=r'Unresolved, '+lmodel_[0],ls='--',marker='o')
	else:
		plt.plot(lz0, lHII0[lu0[3]>0],label=r'Unresolved: radiative, '+lmodel_[1],marker='o')
		plt.plot(lz1, lHII1[lu1[3]>0],label=r'Unresolved: radiative, '+lmodel_[0],ls='--',marker='o')
	if sfdbk==0:
		lHII0_ = epsilon_ff_*lu0[3]*UM/(1e10*mmw()*100*PROTON)
		lHII1_ = epsilon_ff_*lu1[3]*UM/(1e10*mmw()*100*PROTON)
		plt.plot(lz0, lHII0_[lu0[3]>0],label=r'Unresolved: collisional, '+lmodel_[1],lw=1,marker=u'$\u2193$')#r'$\downarrow$')
		plt.plot(lz1, lHII1_[lu1[3]>0],label=r'Unresolved: collisional, '+lmodel_[0],lw=1,ls='--',marker=u'$\u2193$')#'$\downarrow$')
	plt.xlabel(r'$z$')
	plt.ylabel(r'$L_{\mathrm{ff}}\ [\mathrm{erg\ s^{-1}}]$')
	if sca!=0:
		plt.yscale('log')
		#plt.xscale('log')
	plt.legend()
	plt.tight_layout()
	if sca==0:
		plt.savefig(rep0+'L_z_'+str(bins)+'.pdf')
	else:
		plt.savefig(rep0+'logL_z_'+str(bins)+'.pdf')

	"""
	plt.figure()
	plt.plot(out0[0][out0[1]>0],out0[1][out0[1]>0]/out1[1][out0[1]>0],label=r'$L_{\mathrm{WDM}}/L_{\mathrm{CDM}}$')
	plt.plot(out0[0][out0[1]>0],[1.0 for x in range(out0[0][out0[1]>0].shape[0])],color='k',ls='--')
	plt.xlabel(r'$z$')
	plt.ylabel(r'$L_{\mathrm{WDM}}/L_{\mathrm{CDM}}$')
	if sca!=0:
		plt.yscale('log')
		plt.xscale('log')
	plt.tight_layout()
	if sca==0:
		plt.savefig(rep0+'Lrat_z_'+str(bins)+'.pdf')
	else:
		plt.savefig(rep0+'logLrat_z_'+str(bins)+'.pdf')
	#plt.show()
	"""




