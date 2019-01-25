from radio import *

if __name__ == "__main__":
	ncore = 8
	#rep0 = 'halo1_jj_new/'
	rep0 = 'halo1_jj/'
	#rep0 = 'halo3_jj/'
	#rep0 = 'halo1/'
	ldir = ['NL4_zoom_wdm/'+rep0, 'NL4_zoom_cdm/'+rep0]
	#ldir = ['halo1_wdm/','halo1_cdm/']
	nsh = 1.0e-2
	bins = int(sys.argv[1])
	if len(sys.argv)>3:
		low, up = int(sys.argv[2]), int(sys.argv[3])
	else:
		low, up = 1750, 2250

	sn0 = 25
	sn1 = 25
	zmin = 15

	dlu0, dlu1 = [], []
	for sn in range(0,28):
		if sn<=sn0:
			lu0_temp = luminosity_tot(sn,ldir[0],[[low]*3,[up]*3],ncore=ncore,nsh=nsh,zmin=zmin)
			dlu0.append(lu0_temp)
		if sn<=sn1:
			lu1_temp = luminosity_tot(sn,ldir[1],[[low]*3,[up]*3],ncore=ncore,nsh=nsh,zmin=zmin)	
			dlu1.append(lu1_temp)
	lu0 = np.array(dlu0).T
	lu1 = np.array(dlu1).T
	totxt(rep0+'Ltot_z_'+str(bins)+'_'+lmodel[1]+'.txt',lu0,['z', 'L(WDM)', 'delta', 'Msink'],1,0)
	totxt(rep0+'Ltot_z_'+str(bins)+'_'+lmodel[0]+'.txt',lu1,['z', 'L(CDM)', 'delta', 'Msink'],1,0)
