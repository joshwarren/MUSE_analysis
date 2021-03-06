;################################################################
; This is addapted from kinemetry_example.pro for my own use with
; MUSE observations.
;###############################################################

PRO do_work, gal, opt, type
	print, gal, ' ', type
	plot = 0
	; Odd or even moment?
	if stregex(type, '.*vel', /boolean) then even=0 else even=1


	print, 'Have you got the most recent files for '+gal+'?'
	; Select correct file
	file = '/Data/muse/analysis/'+gal+'/'+opt+'/kinemetry/'+type+'.dat'
	; read in field
	rdfloat, file, velbin, er_velbin
	if strcmp(type, 'flux', 4, /FOLD_CASE) then er_velbin = SQRT(velbin)

	; Read in binning
	file = '/Data/muse/analysis/'+gal+'/'+opt+$
		'/setup/voronoi_2d_binning_output.txt'
	rdfloat, file, _,_,bin_num, xbin, ybin, skipline=1

	file = '/Data/muse/analysis/galaxies.txt'
	readcol, file, galaxy_gals,x0,y0, skipline=1,format='A,I,I'
	i_gal = where(galaxy_gals eq gal)
	x0 = float(x0[i_gal[0]])
	y0 = float(y0[i_gal[0]])

	goodpix = where(velbin ne 9999)


	b = uniq(bin_num,sort(bin_num))
	xbin = xbin[b]
	ybin = ybin[b]

	; Center the origin on the center of the galaxy
	xbin = xbin - x0
	ybin = ybin - y0

	xbin = xbin[goodpix]
	ybin = ybin[goodpix]
	velbin = velbin[goodpix]
	er_velbin = er_velbin[goodpix]

	; nor = median(velbin)
	; velbin = velbin / nor
	; er_velbin = er_velbin / nor


	; NB: gas must be the first 3 characters in type
	if strcmp(type, 'gas', 3, /FOLD_CASE) then begin

		; file = '/Data/muse/analysis/'+gal+'/'+opt+'/kinemetry/gas_flux.dat'
		; ; read in field
		; rdfloat, file, flux, er_flux

		; file = '/Data/muse/analysis/'+gal+'/'+opt+$
		; 	'/setup/voronoi_2d_binning_output.txt'
		; rdfloat, file, _,_,bin_num2, xbin2, ybin2, skipline=1
		; goodpix = where(flux ne 9999)

		; b = uniq(bin_num1,sort(bin_num2))
		; xbin2 = xbin2[b]
		; ybin2 = ybin2[b]

		; ; Center the origin on the center of the galaxy
		; xbin2 = xbin2 - x0
		; ybin2 = ybin2 - y0

		; xbin2 = xbin2[goodpix]
		; ybin2 = ybin2[goodpix]
		; flux = flux[goodpix]

		; KINEMETRY, xbin2, ybin2, flux, rad, pa, q, cf, $
		; 	ntrm=6, scale=0.2, $; /FIXCEN, 
		; 	even=1, $
		; 	plot_dir='/Data/muse/analysis/'+gal+'/'+opt+$
		; 	'/kinemetry/kinemetry_gas_flux.jpeg'

		; Get average PA
		KINEMETRY, xbin, ybin, velbin, rad, pa, q, cf, $
			ntrm=6, scale=0.2, $; /FIXCEN, 
			even=even, ERROR=er_velbin, $
			er_pa=er_pa, er_cf=er_cf, er_q=er_q, velkin=velkin, $
			velcirc=velcirc, /bmodel, $;paq=[pa,q], radius=rad/0.2, $
			plot_dir='/Data/muse/analysis/'+gal+'/'+opt+'/kinemetry/kinemetry_'+$
			type+'.jpeg'
	; 	med_pa  = MEDIAN(pa)
	; 	;if even eq 0 then 
	; 	; med_pa = med_pa - 90

	; 	med_q = MEDIAN(q)

	; ; Low harmonic terms only
	; 	KINEMETRY, xbin, ybin, velbin, rad, pa, q, cf, $
	; 		ntrm=6, scale=0.2, $; /FIXCEN, 
	; 		even=even, ERROR=er_velbin, $
	; 		er_pa=er_pa, er_cf=er_cf, er_q=er_q, velkin=velkin, $
	; 		velcirc=velcirc, /bmodel
		
		file = '/Data/muse/analysis/'+gal+'/'+opt+'/kinemetry/kinemetry_'+$
			type+'_2Dmodel.txt'
		forprint2, xbin, ybin, velkin, velcirc, width=200, TEXTOUT=file, /SILENT, $
			comment='   xbin      ybin       velkin    velcirc'

		k1 = SQRT(cf[*,1]^2 + cf[*,2]^2)
		erk1 = (SQRT( (cf[*,1]*er_cf[*,1])^2 + (cf[*,2]*er_cf[*,2])^2 ))/k1

		file = '/Data/muse/analysis/'+gal+'/'+opt+'/kinemetry/kinemetry_'+$
			type+'.txt'
		forprint2, rad, pa, er_pa, 1-q, er_q, k1, erk1, width=200, $
			TEXTOUT = file, /SILENT, comment='  radius(arcsec)      pa(deg)   '+$
			'     err         ellip        err           k1           err'
	endif

	; Include high harmonic terms to reach k5.
	if stregex(type, 'stellar.*', /boolean) then begin
		if stregex(type, '.*_flux', /boolean) then ntrm=10 else ntrm=6

		thisDevice = !D.Name
		Set_Plot, 'Z', /COPY

		Device, Set_Resolution=[1000,1000], Z_Buffer=0
		Erase
		
		KINEMETRY, xbin, ybin, velbin, rad, pa, q, cf, $;x0=x0, y0=y0, $
			ntrm=ntrm, scale=0.2, /FIXCEN, even=even, error=er_velbin, $
			er_pa=er_pa, er_q=er_q, er_cf=er_cf, velkin=velkin, $
			velcirc=velcirc, /bmodel, $;cover=0.05,$
			plot_dir='/Data/muse/analysis/'+gal+'/' + opt + $
			'/kinemetry/kinemetry_'+type+'.jpeg'

		if stregex(type, '.*_vel', /boolean) then begin
			k0 = cf[*,0]
			k1 = SQRT(cf[*,1]^2 + cf[*,2]^2)
			k5 = SQRT(cf[*,5]^2 + cf[*,6]^2)
			k51 = k5/k1
			erk1 = (SQRT( (cf[*,1]*er_cf[*,1])^2 + (cf[*,2]*er_cf[*,2])^2 ))/k1
			erk5 = (SQRT( (cf[*,5]*er_cf[*,5])^2 + (cf[*,6]*er_cf[*,6])^2 ))/k5
			erk51 = ( SQRT( ((k5/k1) * erk1)^2 + erk5^2  ) )/k1 

			file = '/Data/muse/analysis/'+gal+'/' + opt + $
				'/kinemetry/kinemetry_'+type+'.txt'
			forprint2, rad, pa, er_pa, q, er_q, k1, erk1, k51, erk51, width=200, $
				TEXTOUT = file, /SILENT, comment='  radius(arcsec)      pa(deg)  '+$
				'      err         ellip        err           k1           err   '+$
				'       k51         err'
		endif

		if stregex(type, '.*_flux', /boolean) then begin
			b4 = cf[*, 8]
			er_b4 = er_cf[*, 8]
			
			file = '/Data/muse/analysis/'+gal+'/' + opt + $
				'/kinemetry/kinemetry_'+type+'.txt'
			forprint2, rad, pa, er_pa, 1-q, er_q, b4, er_b4, width=200, $
				TEXTOUT = file, /SILENT, $
				comment='  radius(arcsec)      pa(deg)        err         ellip  '+$
				'      err           b4           err'
		endif

		if stregex(type, '.*_vel', /boolean) then begin
			file = '/Data/muse/analysis/'+gal+'/' + opt + $
				'/kinemetry/kinemetry_'+type+'_2Dmodel.txt'
			forprint2, xbin, ybin, velkin, velcirc, width=200, TEXTOUT=file, $
				/SILENT, comment='   xbin      ybin       velkin    velcirc'
		endif

	endif

	; plot=0
	; if keyword_set(plot) then begin
	; 	; plot coeffs.
	; 	thisDevice = !D.Name
	; 	Set_Plot, 'Z', /COPY

	; 	Device, Set_Resolution=[1000,1000], Z_Buffer=0
	; 	Erase
		
	; 	; r = GET_SCREEN_SIZE()
	; 	; window, 1, xsize=r[0]*0.3, ysize=r[1]*0.8
	; 	!p.charsize=3
	; 	!y.style=1
	; 	if type eq 'stellar_vel' then !p.multi=[0,1,4] else !p.multi=[0,1,3]
	; 	!Y.MARGIN=[0,0] ; show plots with shared X axis
	; 	!Y.OMARGIN=[5,3] ; allow for space for the axis labels
	; 	ploterror, rad, pa, er_pa, PSYM=-5, TITLE=gal, xtickformat = '(A1)', $
	; 		YTITLE='!7C!X!N [degrees]', YRANGE=[min(pa),max(pa)]
	; 	ploterror, rad, q, er_q, PSYM=-5, YRANGE=[0,1.1], xtickformat = '(A1)', $
	; 		YTITLE='q'
	; 	ploterror, rad, k1, erk1, PSYM=-5, xtickformat = '(A1)', $
	; 		YTITLE='k1 [km/s]',YRANGE=[0,245]
	; 	if type eq 'stellar_vel' then begin
	; 		ploterror, rad, k51, erk51, PSYM=-5, XTITLE='R [arcsec]', $
	; 		YTITLE='k5/k1', YRANGE=[0,0.13]
	; 	endif
	; 	!P.MULTI=0
	; 	!Y.MARGIN=[4,2] ; back to default values
	; 	!Y.OMARGIN=[0,0]
	; 	!p.charsize=1
	; 	; WAIT, 3000000000


	; 	snapshot = TVRD()
	; 	TVLCT, r, g, b, /Get
	; 	Device, Z_Buffer=1
	; 	Set_Plot, thisDevice

	; 	image24 = BytArr(3, 1000, 1000)
	; 	image24[0,*,*] = r[snapshot]
	; 	image24[1,*,*] = g[snapshot]
	; 	image24[2,*,*] = b[snapshot]

	; 	Write_JPEG, '/Data/muse/analysis/'+gal+'/'+opt+'/kinemetry/kinemetry_'+$
	; 		type+'.jpeg', image24, True=1, Quality=75
	; endif

END


pro use_kinemetry_muse
	do_work, 'ngc1316', 'pop', 'gas_vel'
	; do_work, 'ngc1316', 'pop', 'stellar_vel'
	; gals=['ic1459', 'ic4296', 'ngc1316', 'ngc1399']
	; ; gals=['ngc1316']
	; for i=0,3 do begin
	; 	gal=gals[i]
	; 	do_work, gal, 'kin', 'stellar_flux'
	; 	; do_work, gal, 'kin3', 'stellar_vel'
	; 	; do_work, gal, 'kin3', 'gas_flux'
	; 	; do_work, gal, 'kin3', 'gas_vel'
	; 	; do_work, gal, 'kin3', 'gas_sigma'

	; endfor

END