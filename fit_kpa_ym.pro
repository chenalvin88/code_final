pro test

  rcrit = 0.3
  count = 0
  indxtmp = -1
  samp = mrdfits('/Volumes/SDrive/yenting_pa_alignment/files/dapall-v3_1_1-3.1.0.fits',1)
  readcol,'rg.list',plateifu,format='(a)'
  for ii = 0, n_elements(plateifu)-1 do begin
    plateifu_rg = strtrim(plateifu[ii],2)
    plateifu_mpl11 = strtrim(samp.plateifu,2)
    match,plateifu_rg,plateifu_mpl11,sub1,sub2
    indxtmp = [indxtmp,sub2]
  endfor
  indxtmp = indxtmp[2:n_elements(indxtmp)-1]
  ; indxtmp = where(strtrim(samp.plateifu,2) eq '8154-3702')
  samp = samp[indxtmp]
  nobj = n_elements(samp)

  str = {plateifu: '*', PA_STAR_1RE: 0.D, $
            PA_STAR_ERR_1RE: 0.0D,PA_STAR_VELOFF_1RE: 0.0D,PA_STAR_PXL_1RE: 0.0D,  PA_HA_1RE: 0.0D, $
            PA_HA_ERR_1RE: 0.0D,PA_HA_VELOFF_1RE: 0.0D,PA_HA_PXL_1RE: 0.0D, PA_O3_1RE: 0.0D, $
            PA_O3_ERR_1RE: 0.0D,PA_O3_VELOFF_1RE: 0.0D,PA_O3_PXL_1RE: 0.0D}
 
  ss = make_array(val = str, dim=nobj)

  for kk = 0, nobj-1 do begin
    ss[kk].plateifu = samp[kk].plateifu
    print,'the',kk,'th object',samp[kk].plate,' ',samp[kk].ifudesign
    dapname = '/Volumes/SDrive/yenting_pa_alignment/MaNGA/SPX-MILESHC-MASTARSSP/'+'manga-'+strtrim(string(samp[kk].plate),2)+'-'+strtrim(samp[kk].ifudesign,2)+'-MAPS-SPX-MILESHC-MASTARSSP.fits.gz'
    
    resulttmp = file_search(dapname)
    if strlen(resulttmp) le 3 then continue
        
    nfiber = STRSPLIT(samp[kk].ifudesign, '0',/EXTRACT)
    if strtrim(nfiber[0],2) eq '127' then crit_pixel = 50
    if strtrim(nfiber[0],2) eq '91' then crit_pixel = 50
    if strtrim(nfiber[0],2) eq '61' then crit_pixel = 50
    if strtrim(nfiber[0],2) eq '37' then crit_pixel = 30
    if strtrim(nfiber[0],2) eq '19' then crit_pixel = 20
    crit_pixel = 0.5

    ;;;;;;;read DAP results;;;;;;;;;
    vstar = mrdfits(dapname, 15)
    vstarivar = mrdfits(dapname,16)
    vstarmask = mrdfits(dapname,17)

    vgas = mrdfits(dapname,37)
    vgasivar = mrdfits(dapname,38)
    vgasmask = mrdfits(dapname,39)

    vha = vgas[*,*,23]
    vhaivar = vgasivar[*,*,23]
    vhamask = vgasmask[*,*,23]

    vo3 = vgas[*,*,16]
    vo3ivar = vgasivar[*,*,16]
    vo3mask = vgasmask[*,*,16]

    vstarbak = vstar
    vhabak = vha
    vo3bak = vo3

    fgas = mrdfits(dapname,30)
    fgasivar = mrdfits(dapname,31)
    fgasmask = mrdfits(dapname,32)

    iok = where(fgasivar GT 0)
    fgaserr = 0 * fgasivar
    if (iok[0] NE -1) then fgaserr[iok] = 1 / sqrt(fgasivar[iok])

    fha = fgas[*,*,23]
    fhaerr = fgaserr[*,*,23]

    fo3 = fgas[*,*,16]
    fo3err = fgaserr[*,*,16]

    si = size(vstar)
    xpos = dblarr(si[1], si[2])
    ypos = xpos
    nx = si[1]
    ny = si[2]
    nc = fix(nx/2.)

    for xx = 0, nx-1 do begin
      for yy = 0, ny-1 do begin
        xpos[xx,yy] = (xx - nx/2.0 + 0.5)/2.
        ypos[xx,yy] = (yy - ny/2.0 + 0.5)/2.
      endfor
    endfor

    x = xpos
    y = ypos

    r = sqrt(x^2 + y^2)*2
    r = minmax(sqrt(x^2+y^2)*(2*(x ge 0)-1))

    del = (-90+samp[kk].NSA_ELPETRO_PHI)/180.*!pi
    ba = samp[kk].NSA_ELPETRO_BA
    p50 = samp[kk].NSA_elpetro_th50_r
    ; print, 'del',del,'ba',ba,'p50',p50
    xpos2 = xpos*cos(del)+ypos*sin(del)
    ypos2 = -xpos*sin(del)+ypos*cos(del)
    r0 = sqrt(xpos2^2+(ypos2/ba)^2)
    r1 = r0/p50
    indstar = where(vstarmask ne 1 and vstarivar gt 0 and r1 lt rcrit,ctstar)
    ; print, 'vstarmask',vstarmask
    ; print, 'vstarivar',vstarivar
    ; print, 'r1',r1
    ; print, 'indstar',size(indstar),indstar
    indha = where(vhamask ne 1 and vhaivar gt 0 and fha/fhaerr gt 3 and r1 lt rcrit,ctha)
    indo3 = where(vo3mask ne 1 and vo3ivar gt 0 and r1 lt rcrit and fo3/fo3err gt 3,cto3)

    v0 = vstar[nc,nc]
    ; print, v0
    vstar = vstar - v0 ; IMPORTANT: Subtract v0 before the fit
    xstar = x[indstar]
    ystar = y[indstar]
    vstar = vstar[indstar]

    if ctha gt crit_pixel then begin
      v0 = vha[nc,nc]
      vha = vha - v0 ; IMPORTANT: Subtract v0 before the fit
      xha = x[indha]
      yha = y[indha]
      vha = vha[indha]
    endif

    if cto3 gt crit_pixel then begin
      v0 = vo3[nc,nc]
      vo3 = vo3 - v0 ; IMPORTANT: Subtract v0 before the fit
      xo3 = x[indo3]
      yo3 = y[indo3]
      vo3 = vo3[indo3]
    endif
    ; print, 'xstar',xstar,'ystar',ystar,'vstar',vstar
    if ctstar gt crit_pixel then fit_kinematic_pa, xstar, ystar, vstar, pastar, pastarErr, VelSyststar, NSTEPS=30
    ; print, 'star',pastar, pastarErr, VelSyststar,ctstar
    ; print, 'ctha',ctha
    if ctha gt crit_pixel then fit_kinematic_pa, xha, yha, vha, paha, pahaErr, VelSystha, NSTEPS=30
    ; print, 'ha',paha,pahaErr,VelSystha,ctha
    ; print, 'cto3',cto3
    if cto3 gt crit_pixel then fit_kinematic_pa, xo3, yo3, vo3, pao3, pao3Err, VelSysto3, NSTEPS=30
    print, 'o3',pao3, pao3Err, VelSysto3,cto3

    if ctstar gt crit_pixel then begin
      ss[kk].PA_STAR_1RE = pastar
      ss[kk].PA_STAR_ERR_1RE = pastarErr
      ss[kk].PA_STAR_VELOFF_1RE = VelSyststar
      ss[kk].PA_STAR_PXL_1RE = ctstar
    endif

    if ctha gt crit_pixel then begin
      ss[kk].PA_HA_1RE = paha
      ss[kk].PA_HA_ERR_1RE = pahaErr
      ss[kk].PA_HA_VELOFF_1RE = VelSystha
      ss[kk].PA_HA_PXL_1RE = ctha
    endif

    if cto3 gt crit_pixel then begin
      ss[kk].PA_O3_1RE = pao3
      ss[kk].PA_O3_ERR_1RE = pao3Err
      ss[kk].PA_O3_VELOFF_1RE = VelSysto3
      ss[kk].PA_O3_PXL_1RE = cto3
    endif
    write_csv,'result_ym.csv',ss
  endfor

  ; mwrfits,ss,'MPL11_rg_PA_1Re_for_Yenting.fits',/create

end
