;################################################################
; This is an example routine which calls KINEMETRY 
; to analyse SAURON velocity map of NGC2974 as presented in
; (Emsellem et al. 2004 MNRAS, 352, 271). It makes a plot
; of kinemetric coefficients as in the first column of Fig.7 
; in Krajnovic et al. (2006).
;
; This routine uses RDFLOAT.PRO and PLOTERROR.PRO avaiable from
; IDL Astronomer Unser's Library.
; 
; Davor Krajnovic, Oxford, 07.12.2005.
;###############################################################

PRO kinemetry_example

    ;
    ; read in all data
    ;
    file = './kinemetry/NGC2974_SAURON_kinematics.dat'
    rdfloat, file, num, xbin, ybin, velbin, er_velbin,  SKIPLIN=1
    readcol,'rg.list',plateifu,format='(a)'
    radlist = make_array(35,n_elements(plateifu),/string,value=0)
    palist = make_array(35,n_elements(plateifu),/string,value=0)
    paerrlist = make_array(35,n_elements(plateifu),/string,value=0)
    qlist = make_array(35,n_elements(plateifu),/string,value=0)
    qerrlist = make_array(35,n_elements(plateifu),/string,value=0)
    k1list = make_array(35,n_elements(plateifu),/string,value=0)
    k1errlist = make_array(35,n_elements(plateifu),/string,value=0)
    k51list = make_array(35,n_elements(plateifu),/string,value=0)
    k51errlist = make_array(35,n_elements(plateifu),/string,value=0)

    ; for i = 0, 0 do begin
    for i = 0, n_elements(plateifu)-1 do begin
        print, plateifu[i]
        dapname = '/Volumes/SDrive/yenting_pa_alignment/MaNGA/VOR10-MILESHC-MASTARSSP/manga-'+plateifu[i]+'-MAPS-VOR10-MILESHC-MASTARSSP.fits.gz'
        coo = mrdfits(dapname, 1)
        X = -coo[*,*,0]
        Y = coo[*,*,1]
        ; print, x,y
        vstar = mrdfits(dapname, 15)
        vstarivar = mrdfits(dapname,16)
        vstarmask = mrdfits(dapname,17)
        indstar = where(vstarmask ne 1,ctstar)
        xbin = X[indstar]
        ybin = Y[indstar]
        velbin = vstar[indstar]
        er_velbin = 1/sqrt(vstarivar)
        maxind = floor(max([xbin,ybin]))

        ;
        ; kinemetry on velocity map
        ;
        t=systime(1)
        KINEMETRY, xbin, ybin, velbin, rad, pa, q, cf, ntrm=6, scale=0.5, $
            ERROR=er_velbin, name='NGC2974',er_cf=er_cf, er_pa=er_pa, $
            er_q=er_q, /plot, /verbose
        ; print, systime(1) -t, 'seconds'

        ;
        ; kinemetry parameters as defined in Krajnovic et al. (2006)
        ; 
        k0 = cf[*,0]
        k1 = SQRT(cf[*,1]^2 + cf[*,2]^2)
        k5 = SQRT(cf[*,5]^2 + cf[*,6]^2)
        k51 = k5/k1
        erk1 = (SQRT( (cf[*,1]*er_cf[*,1])^2 + (cf[*,2]*er_cf[*,2])^2 ))/k1
        erk5 = (SQRT( (cf[*,5]*er_cf[*,5])^2 + (cf[*,6]*er_cf[*,6])^2 ))/k5
        erk51 = ( SQRT( ((k5/k1) * erk1)^2 + erk5^2  ) )/k1 

        ;
        ; plot coeffs.
        ;
        ; r = GET_SCREEN_SIZE()
        ; window, 1, xsize=r[0]*0.3, ysize=r[1]*0.8
        ; !p.charsize=3
        ; !y.style=1
        ; !p.multi=[0,1,4]
        ; !Y.MARGIN=[0,0] ; show plots with shared X axis
        ; !Y.OMARGIN=[5,3] ; allow for space for the axis labels
        ; ploterror, rad, pa, er_pa, PSYM=-5, TITLE='NGC2974', xtickformat = '(A1)', YTITLE='!7C!X!N [degrees]', YRANGE=[30,70], /Current
        ; ploterror, rad, q, er_q, PSYM=-5, YRANGE=[0,1.1], xtickformat = '(A1)', YTITLE='q'
        ; ploterror, rad, k1, erk1, PSYM=-5, xtickformat = '(A1)', YTITLE='k1 [km/s]',YRANGE=[0,245]
        ; ploterror, rad, k51, erk51, PSYM=-5, XTITLE='R [arcsec]', YTITLE='k5/k1', YRANGE=[0,0.13]
        ; !P.MULTI=0
        ; !Y.MARGIN=[4,2] ; back to default values
        ; !Y.OMARGIN=[0,0]
        ; !p.charsize=1

        n=n_elements(rad)
        if n le maxind then maxind=n-1
        !P.Region=[0,0,1,1]
        Plot = errorplot(rad[0:maxind],pa[0:maxind],er_pa[0:maxind],$
            POSITION=[.15,.70,0.95,.90], TITLE=plateifu[i],xstyle=1, xmajor=0, YTITLE='PA [degrees]',ystyle=2)
        Plot = errorplot(rad[0:maxind],q[0:maxind], er_q[0:maxind],$
            /CURRENT, POSITION=[.15,.50,0.95,.70],xstyle=1, xmajor=0, YTITLE='q',ystyle=2)
        Plot = errorplot(rad[0:maxind],k1[0:maxind], erk1[0:maxind],$
            /CURRENT, POSITION=[.15,.30,0.95,.50],xstyle=1, xmajor=0, YTITLE='k1 [km/s]',ystyle=2)
        Plot = errorplot(rad[0:maxind],k51[0:maxind], erk51[0:maxind],$
            /CURRENT, xtickvalues = rad[0:maxind],POSITION=[.15,.10,0.95,.30],xstyle=1,YTITLE='k5/k1',yrange=[0,0.2])
        aPlot = errorplot(rad[0:maxind],k51[0:maxind], erk51[0:maxind], /OVERPLOT, xmajor=0,yrange=[0,0.2])
        aPlot.Save, '/Volumes/SDrive/yenting_pa_alignment/results/kinemetry/'+plateifu[i]+'.png'
        aPlot.close
        radlist[0:n,i]=[plateifu[i],string(rad)]
        palist[0:n,i]=[plateifu[i],string(pa)]
        paerrlist[0:n,i]=[plateifu[i],string(er_pa)]
        qlist[0:n,i]=[plateifu[i],string(q)]
        qerrlist[0:n,i]=[plateifu[i],string(er_q)]
        k1list[0:n,i]=[plateifu[i],string(k1)]
        k1errlist[0:n,i]=[plateifu[i],string(erk1)]
        k51list[0:n,i]=[plateifu[i],string(k51)]
        k51errlist[0:n,i]=[plateifu[i],string(erk51)]
        ; ['pa',string(plateifu[i]),string(pa)]
        write_csv,'/Volumes/SDrive/yenting_pa_alignment/results/kinemetry/result_rad.csv',radlist
        write_csv,'/Volumes/SDrive/yenting_pa_alignment/results/kinemetry/result_pa.csv',palist
        write_csv,'/Volumes/SDrive/yenting_pa_alignment/results/kinemetry/result_paerr.csv',paerrlist
        write_csv,'/Volumes/SDrive/yenting_pa_alignment/results/kinemetry/result_q.csv',qlist
        write_csv,'/Volumes/SDrive/yenting_pa_alignment/results/kinemetry/result_qerr.csv',qerrlist
        write_csv,'/Volumes/SDrive/yenting_pa_alignment/results/kinemetry/result_k1.csv',k1list
        write_csv,'/Volumes/SDrive/yenting_pa_alignment/results/kinemetry/result_k1err.csv',k1errlist
        write_csv,'/Volumes/SDrive/yenting_pa_alignment/results/kinemetry/result_k51.csv',k51list
        write_csv,'/Volumes/SDrive/yenting_pa_alignment/results/kinemetry/result_k51err.csv',k51errlist

    endfor
END
