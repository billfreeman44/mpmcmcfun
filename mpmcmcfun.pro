; NAME:
;   mpmcmcfun
;
; PURPOSE:
;   This program is used to fit any function using Markov Chain Monte Carlo (MCMC)
;   The idea behind it is that its supposed to be like "mpfitfun" but with MCMC.
;   Hence the name...
;
;
; 
; The program requires the Coyote Library (http://www.idlcoyote.com/documents/programs.php)
; to be installed on your machine.  Also, the two programs likelih_mpmcmcfun and 
; calc_mom_mpmcmcfun
;
;
;RETURN VALUE:
; The best fit values are returned.  The number 0 may be returned in some cases if failures occured.
; Checking status_MCMC is recommended.
;
;    
;PARAMETERS:
;All parameters are required
;    MYFUNCT: String of the function you wish to fit.  For instance if you 
;     wanted to fit a gaussian this would just be "gaussian"
;    
;    X: X array of data you wish to fit
;    
;    Y: Y array of data you wish to fit
;    
;    ERR: One sigma errors of your Y data
;    
;    P: Initial guess of your parameters.  I recommend running the real
;      mpfitfun on your initial data then plugging in the output to this function
;    
;    P_SIGMA: Sigma of your guesses.  This is a very important input.  If you don't know
;      what this is, read up on the Metropolis-Hastings algorithm.  All values must be 
;      positive and non-zero
;    
;    n_iter: Number of iterations.
;       
;KEYWORDS
;
;  parinfo: Just like parinfo in mpfitfun, although this parinfo CANNOT TAKE EVERY PARAMETER
;     THAT THE NORMAL PARINFO IN MPFIT CAN.  It can only do these three:
;     {fixed:0, limited:[1,1], limits:[0.0D,0.0d], tied:''}
;     
;  PERROR: One sigma errors in the fit.  This is usually the one sigma fit to the gaussian histogram.
;     if there is no good fit, then 
;     
;  calculated_params: An array of strings of things to calculate using P[#] for each parameter.  For
;     instance: ["p[1]*p[2]","p[4]*p[2]/P[6]"] would work.  It must be an array.
;     
;  titles: An array of strings to be the plot titles.  THIS MUST INCLUDE THE CALCULATED PARAMS AS WELL. See
;     the example program.
;     
;  makeplots: Flag for weather or not to make plots.  Uses cgplot.  You should have run ps_start before.
;     see the example as well.
;     
;  status_MCMC: output.  1 if the fit went ok, 0 if it failed.  I recommend you check this.
;  
;  gauss_fit_worked: output array in parallel with your guesses and calculated parameters.  1  if the gaussian
;     fit worked and 0 if it didin't.
;     
;  makecontours: Flag for weather or not to make contours.  Uses cgcontour.  You should have run ps_start before.
;     see the example as well.  Makes contours for each parameter.
;  
;  silent or quiet: set either to not output any text
;  
;  n_bins: number of bins used in the histogram (int) (default 50, optional)
;
;  cup_pars: parameters calculated by "filling the cup"
;    0 best fit 
;    1 +1sigma
;    2 -1sigma
;    3 +2sigma
;    4 -2sigma
;    5 +3sigma
;    6 -3sigma
;    
;  plotchain: keyword. set to plot chain.
;          
;EXAMPLES
;    See mpmcmcfun_test.pro
;    
;       
;        
;AUTHOR
;       Bill Freeman:
;        Billfreeman44@yahoo.com
;
;
;
;
;
;Copyright (c) 2014 William R. Freeman
;
;Permission is hereby granted, free of charge, to any person 
;obtaining a copy of this software and associated documentation
;files (the "Software"), to deal in the Software without restriction, 
;including without limitation the rights to use, copy, modify, merge,
;publish, distribute, sublicense, and/or sell copies of the 
;Software, and to permit persons to whom the Software is furnished 
;to do so, subject to the following conditions:
;
;The above copyright notice and this permission 
;notice shall be included in all copies or substantial 
;portions of the Software.
;
;THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY 
;OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
; TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
; PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
; THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
; DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
; CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
; IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
; DEALINGS IN THE SOFTWARE.
;
;
;
;;git commit -a -m ''
;result=mpfitfun('mp_sngle_gauss_n2',wavel[index],data[index],err[index],p,parinfo=pi,/QUIET,status=status,DOF=DOF,PERROR=PERROR1)
;parinfo should have these 3: {fixed:0, limited:[1,1], limits:[0.D,1.2d], tied:''}
;calculated_params could be like ["p[1]*p[2]","p[4]*p[2]/P[6]"]
function mpmcmcfun,MYFUNCT, X, Y, ERR, P,P_SIGMA,n_iter ,parinfo=parinfo,$
  PERROR=PERROR,calculated_params=calculated_params,titles=titles,makeplots=makeplots,status_MCMC=status_MCMC,$
  gauss_fit_worked=gauss_fit_worked,makecontours=makecontours,silent=silent,quiet=quiet,$
  format_str=format_str,n_bins=n_bins,cup_pars=cup_pars,plotchain=plotchain,$
  debug=debug
if n_params() ne 7 then message,'ALL PARAMETERS ARE REQUIRED. syntax: result=mpfitfun_mcmc(MYFUNCT, X, Y, ERR,'+$
  ' P,P_SIGMA,n_iter ,parinfo,PERROR=PERROR,calculated_params=calculated_params,titles=titles'+$
  ',makeplots=makeplots,status_MCMC=status_MCMC,$'+$
  'gauss_fit_worked=gauss_fit_worked,makecontours=makecontours)'
  

;SET UP TIME COUNTER AND DEFAULT KEYWORDS.
TIMEvar=systime(/seconds)
status_MCMC=1
if keyword_set(quiet) then silent=1
if ~keyword_set(calculated_params) then begin
  calculated_params=[]
  calculate_pars=0
  endif else calculate_pars=1
if ~keyword_set(titles) then titles=replicate(' ',n_elements(P)+n_elements(calculated_params))
if ~keyword_set(parinfo) then parinfo=replicate({fixed:0, limited:[0,0], limits:[0.D,0.d], tied:''},n_elements(p))
if ~keyword_set(n_bins) then n_bins=50


;ensure double precision
xin=double(X)
yin=double(Y)
yinerr=double(ERR)
seed=1234d
guess=P
guess_old=guess
guess_arr=replicate(0.d,n_elements(guess),n_iter)
accepted_arr=replicate(-1l,n_elements(guess),n_iter)


;initial guess likelihood
t=execute('yout='+MYFUNCT+'(xin,guess)')
if t eq 0 then message,'bad MYFUNCT'
l=likelih_mpmcmcfun(xin,yin,yinerr,yout)
l_old=10^double(l)
if l_old eq 0 then begin
  print,'initial guess likleihood is 0, check initial parameters.'
  status_MCMC=0
  print,l

;  close,lun
;  free_lun,lun
  return,p
  endif

;calc which parameters are free. guess_small has only free parameters where guess has all of them.
;;{fixed:0, limited:[1,1], limits:[0.D,0.5d], tied:''}
guess_small_indicies=where(parinfo.fixed eq 0 and (parinfo.tied eq '' or parinfo.tied eq ' '))
guess_small=guess[guess_small_indicies]
guess_small_old=guess_small
fixed_pars=where(parinfo.tied ne '' and parinfo.tied ne ' ',/null)

;check for 0 sigma
for i=0,n_elements(guess_small_indicies)-1 do $
  if P_SIGMA[guess_small_indicies[i]] le 0 then begin
    status_MCMC=0
    if ~keyword_set(silent) then print,'WARNING: negative or zero sigma.'
    return,P
    endif

;counter setup
n_better=0l
n_worse=0l
np=float(n_elements(guess_small))

;======================================================================
;ENTER MAIN LOOP OF MCMC.
;======================================================================
seed=1234d
for i=0l,n_iter-1 do begin
  ;change one random parameter by random amount 
  change_index=fix(randomu(seed)*np)
  guess_small[change_index]=guess_small_old[change_index]+randomn(seed)*P_SIGMA[guess_small_indicies[change_index]]
;  di=guess[guess_small_indicies[change_index]]
  ;adjust main guess by small guess.
  guess[guess_small_indicies]=guess_small
  p=guess
  for j=0,n_elements(fixed_pars)-1 do t=execute('guess[fixed_pars[j]]='+parinfo[fixed_pars[j]].tied)
  
  ;check things. limited, limits
  limit_exceeded=0
  for k=0,n_elements(guess)-1 do begin
    if (parinfo[k].tied eq '' or parinfo[k].tied eq ' ') AND PARINFO[K].fixed ne 1 then begin 
      if parinfo[k].limited[0] eq 1 then if guess[k] lt parinfo[k].limits[0] then limit_exceeded=1
      if parinfo[k].limited[1] eq 1 then if guess[k] gt parinfo[k].limits[1] then limit_exceeded=1
      endif;fixed parinfo
    endfor   

  if limit_exceeded eq 0 then begin
    YOUT=call_function(myfunct,xin,guess)
    l=likelih_mpmcmcfun(xin,yin,yinerr,yout)
    l=2.7182818^double(l)
    if keyword_set(debug) then print,'new, old l,ratio: ',l_old,l,l/l_old
    endif else l=0.0
  
;  if change_index eq 0 then print, l / l_old,guess_small[change_index],guess_small[change_index]-di;,change_index
  if l / l_old gt randomu(seed) then begin
    ;fit is better
    guess_old=guess
    guess_small_old=guess[guess_small_indicies]
    l_old=l
    n_better=n_better+1l  
;    print,i,n_better,n_worse+n_better
    guess_arr[*,n_better-1]=guess ; add to chain
    accepted_arr[guess_small_indicies[change_index],i]=1
    endif else begin
      ;fit is worse
      guess=guess_old
      guess_small=guess_small_old
      n_worse=n_worse+1l
;      print,i,n_worse,n_worse+n_better
      accepted_arr[guess_small_indicies[change_index],i]=0
      endelse  
  endfor ; n_iter


;print how many of each parameter are accepted.
if ~keyword_set(silent) then begin
  print,'n better ',n_better,double(n_better)/double(n_iter)
  print,'n worse  ', n_worse,double(n_worse)/double(n_iter)
  tote=0l
  for i=0,n_elements(guess)-1 do begin
    index=where(accepted_arr[i,*] ne -1)
    ind_good=where(accepted_arr[i,index] eq 1,ctg)
    ind_bad=where(accepted_arr[i,index] eq 0,ctb)
    tote=tote+ctg+ctb
    if ctb ne 0 then print,i,titles[i],ctg,ctb,double(ctg)/double(ctb+ctg),P_SIGMA[i] 
    endfor
  endif


;check if low number of successful/failed jumps
if n_better lt 0.01*n_iter or n_worse lt 0.01*n_iter then begin
  help
  print,'not enough better/worse soultions.'
  print,'n better ',n_better,double(n_better)/double(n_iter)
  print,'n worse  ', n_worse,double(n_worse)/double(n_iter)
  tote=0l
  for i=0,n_elements(guess)-1 do begin
    index=where(accepted_arr[i,*] ne -1)
    ind_good=where(accepted_arr[i,index] eq 1,ctg)
    ind_bad=where(accepted_arr[i,index] eq 0,ctb)
    tote=tote+ctg+ctb
    if ctb ne 0 then print,i,titles[i],ctg,ctb,double(ctg)/double(ctb+ctg),P_SIGMA[i] 
    endfor
  status_MCMC=0
;  stop
  return,p
  endif

;======================================================================
;ENTER MAIN PLOTTING AND ANALYSIS SECTION.
;======================================================================
return_value=[]
return_value_errs=[]
cup_pars=replicate(-1.0d,n_elements(guess)+n_elements(calculated_params),7)
gauss_fit_worked=[]

;create arrays of the chain!
for i=0,n_elements(guess)-1 do t=execute("CH"+SSI(I)+"=guess_arr["+ssi(i)+",0:n_better-1]")

for i=0,n_elements(guess)-1 do begin
  ;check if data was not fit and fixed to constant.
  skipguess=0
  t=execute("if min(CH"+SSI(I)+") - max(CH"+ssi(i)+") eq 0 then skipguess=1")
  if skipguess eq 1 then goto,skipguess_location
  
  ;calc histogram
  if keyword_set(makeplots) then begin
  !p.multi=[0,1,1]
  
  ;fix infinite values of ch
;  if titles[i] eq 'OIII/Hbeta broad' then begin
;    t=execute("CH"+SSI(I)"=  ")
;    endif
    
    if (parinfo[i].tied eq '' or parinfo[i].tied eq ' ') AND PARINFO[i].fixed ne 1 then $
      t=execute("cghistoplot,CH"+SSI(I)+$
      ",nbins="+ssi(n_bins)+",title=titles[i]+' '+ssf(P_SIGMA[i]),histdata=histdata,locations=locations") $
      else t=execute("histdata=histogram(CH"+SSI(I)+",nbins="+ssi(n_bins)+",locations=locations)")
    endif else t=execute("histdata=histogram(CH"+SSI(I)+",nbins="+ssi(n_bins)+",locations=locations)")
  
  ;calculate initial guess and initial conditions.
  bestval=calc_mom_mpmcmcfun(locations,histdata,1)
  sigval=calc_mom_mpmcmcfun(locations,histdata,2)
  pi2 =[{fixed:0, limited:[1,1], limits:[0.D,max(histdata)*1.5]},$ ;peak value
    {fixed:0, limited:[0,0], limits:[0.D,0.d]},$ ;peak centroid
    {fixed:0, limited:[0,0], limits:[0.D,0.D]}] ;sigma
  estimates=[max(histdata),bestval,sigval]
  
  ;do the fit
  dummy=MPFITPEAK(locations,histdata,$
    coeff,nterms=3,/gaussian,$
    estimates=estimates,parinfo=pi2,status=status,chisq=chisq)
  
  
  if (parinfo[i].tied eq '' or parinfo[i].tied eq ' ') AND PARINFO[i].fixed ne 1 then begin 
    if parinfo[i].limited[0] eq 1 then if coeff[1] lt parinfo[i].limits[0] then coeff[1]=parinfo[i].limits[0]
    if parinfo[i].limited[1] eq 1 then if coeff[1] gt parinfo[i].limits[1] then coeff[1]=parinfo[i].limits[1]
    endif;fixed parinfo
  
  ;cup_pars[par number,max/sigma]
  ;calculate "cup" parameters
  cup_out=mpmcmc_cup(locations,histdata)
  for ii=0,6 do cup_pars[i,ii]=cup_out[ii]
  
  if keyword_set(makeplots) and $
     (parinfo[i].tied eq '' or parinfo[i].tied eq ' ') AND $
     PARINFO[i].fixed ne 1 then $
     for ii=0,6 do cgplot,[cup_out[ii],cup_out[ii]],minmax(histdata),/overplot
  
  
  
  

  skipguess_location:
  if skipguess eq 1 then begin
    t=execute("return_value=[return_value,min(CH"+SSI(I)+")]")
    return_value_errs=[return_value_errs,0.0]
    gauss_fit_worked=[gauss_fit_worked,0]
    endif else begin 
      if status lt 1 then begin
        IF KEYWORD_SET(makeplots)  and $
      (parinfo[i].tied eq '' or parinfo[i].tied eq ' ') AND $
      PARINFO[i].fixed ne 1 THEN BEGIN
          cgplot,locations,gaussian(locations,estimates),/overplot,linestyle=2
          cgtext,0.2,0.9,ssf(bestval),/normal
          cgtext,0.2,0.85,ssf(sigval),/normal
          ENDIF
        return_value=[return_value,bestval]
        return_value_errs=[return_value_errs,sigval]
        gauss_fit_worked=[gauss_fit_worked,0]
        endif else begin
          IF KEYWORD_SET(makeplots)  and $
      (parinfo[i].tied eq '' or parinfo[i].tied eq ' ') AND $
      PARINFO[i].fixed ne 1 THEN BEGIN
            cgplot,locations,dummy,/overplot
            cgplot,locations,gaussian(locations,estimates),/overplot,linestyle=2
            cgtext,0.2,0.9,ssf(coeff[1]),/normal
            cgtext,0.2,0.85,ssf(coeff[2]),/normal
            cgtext,0.2,0.75,ssf(bestval),/normal
            cgtext,0.2,0.70,ssf(sigval),/normal
            if parinfo[i].limited[0] eq 1 then cgtext,0.8,0.90,ssf(parinfo[i].limits[0]),/normal
            if parinfo[i].limited[1] eq 1 then cgtext,0.8,0.85,ssf(parinfo[i].limits[1]),/normal
            ENDIF
          return_value=[return_value,coeff[1]]
          return_value_errs=[return_value_errs,coeff[2]]
          gauss_fit_worked=[gauss_fit_worked,1]
        endelse ;not status ne 1
      
      if keyword_set(plotchain) and $
      (parinfo[i].tied eq '' or parinfo[i].tied eq ' ') AND $
      PARINFO[i].fixed ne 1 then begin
        if n_iter gt 1000002 then $
        t=execute("cgplot,CH"+SSI(I)+"[0:100000],/ynozero") $
          else t=execute("cgplot,CH"+SSI(I)+",/ynozero")
        t=execute("cgplot,[0,1e7],[cup_out[0],cup_out[0]],/overplot")
        t=execute("cgplot,[0,1e7],[cup_out[1],cup_out[1]],/overplot")
        t=execute("cgplot,[0,1e7],[cup_out[2],cup_out[2]],/overplot")
        endif
      
      ;;CONTOURS.
      ;plot all contours.
      if keyword_set(makecontours) then begin
        for j=i+1,n_elements(guess)-1 do begin
          skipcontour=0
          t=execute("if min(CH"+SSI(I)+") - max(CH"+ssi(i)+") eq 0 then skipcontour=1")
          t=execute("if min(CH"+SSI(J)+") - max(CH"+ssi(J)+") eq 0 then skipcontour=1")
          
          if parinfo[j].tied ne '' and parinfo[j].tied ne ' ' then goto,skipcontour_location
          if parinfo[j].fixed eq 1 then goto,skipcontour_location
          if parinfo[i].tied ne '' and parinfo[i].tied ne ' ' then goto,skipcontour_location
          if parinfo[i].fixed eq 1 then goto,skipcontour_location
          
          
;          PRINT,SKIPCONTOUR
          if skipcontour eq 1 then goto,skipcontour_location
          t=execute("xr=minmax(CH"+SSI(I)+")")
          t=execute("yr=minmax(CH"+SSI(j)+")")
;          t=execute("cgplot,ch"+SSI(I)+"[0:5000],ch"+SSI(J)+"[0:5000],xtitle=titles[i],ytitle=titles[j],psym=-6,/ynozero,/ys,/xs,xr=xr,yr=yr,symsize=0.1")
;          if titles[j] eq 'broad ha amp (FREE)' then BEGIN
;            t=execute("cgplot,ch"+SSI(I)+"[5000:15000],ch"+SSI(J)+"[5000:15000],xtitle=titles[i],ytitle=titles[j],psym=6,/ynozero,/ys,/xs,xr=xr,yr=yr,symsize=0.1")
;            t=execute("cgplot,[0,1],[1,0]+0.01,/OVERPLOT")
;            t=execute("cgplot,[0,1],[1,0]+0.00,/OVERPLOT")
;            t=execute("cgplot,[0,1],[1,0]-0.01,/OVERPLOT")
;            t=execute("cgplot,ch"+SSI(I)+"[15000:25000],ch"+SSI(J)+"[15000:25000],xtitle=titles[i],ytitle=titles[j],psym=6,/ynozero,/ys,/xs,xr=xr,yr=yr,symsize=0.1")
;            t=execute("cgplot,[0,1],[1,0]+0.01,/OVERPLOT")
;            t=execute("cgplot,[0,1],[1,0]+0.00,/OVERPLOT")
;            t=execute("cgplot,[0,1],[1,0]-0.01,/OVERPLOT")
;            t=execute("cgplot,ch"+SSI(I)+"[55000:85000],ch"+SSI(J)+"[55000:85000],xtitle=titles[i],ytitle=titles[j],psym=6,/ynozero,/ys,/xs,xr=xr,yr=yr,symsize=0.1")
;            t=execute("cgplot,[0,1],[1,0]+0.01,/OVERPLOT")
;            t=execute("cgplot,[0,1],[1,0]+0.00,/OVERPLOT")
;            t=execute("cgplot,[0,1],[1,0]-0.01,/OVERPLOT")
            
;            ENDIF
          t=execute("binx=(max(ch"+SSI(I)+")-min(ch"+SSI(I)+"))/15.0")
          t=execute("biny=(max(ch"+SSI(j)+")-min(ch"+SSI(j)+"))/15.0")
          t=execute("density=hist_2d(ch"+SSI(I)+",ch"+SSI(j)+",bin1=binx,bin2=biny,min1=min(ch"+SSI(I)+"),max1=max(ch"+SSI(I)+$
            "),min2=min(ch"+SSI(j)+"),max2=max(ch"+SSI(j)+"))")
          
          ;calcluate bin locations
          t=execute("x_locations=[min (ch"+SSI(I)+")]")
          while n_elements(x_locations) lt n_elements(density[*,0]) do x_locations=[x_locations,x_locations[-1]+binx]
          t=execute("y_locations=[min (ch"+SSI(J)+")]")
          while n_elements(y_locations) lt n_elements(density[0,*]) do y_locations=[y_locations,y_locations[-1]+biny]
          x_locations=x_locations+binx/2
          y_locations=y_locations+biny/2
          contour,density,x_locations,y_locations,/ys,/xs,xr=xr,yr=yr,xtitle=titles[i],ytitle=titles[j]
;          t=execute("cgplot,CH"+SSI(I)+",CH"+SSI(j)+",xtitle=titles[i],ytitle=titles[j]")
            skipcontour_location:
            
          endfor;contour j
        endif;contour keyword.  
    ENDELSE ; not skipguess
  endfor;all guesses



;analyze calculated parameters (just like the parameters)
for i=0,n_elements(calculated_params)-1 do BEGIN
  skipcalcp=0
  calculate_par_str=STRLOWCASE(calculated_params[i])
  p_position=strpos(calculate_par_str,'p')
  while p_position ne -1 do begin
    calculate_par_str=strmid(calculate_par_str,0,p_position)+"ch"+strmid(calculate_par_str,p_position+1)
    p_position=strpos(calculate_par_str,'p')
    endwhile
  br_position=strpos(calculate_par_str,'[')
  while br_position ne -1 do begin
    calculate_par_str=strmid(calculate_par_str,0,br_position)+strmid(calculate_par_str,br_position+1)
    br_position=strpos(calculate_par_str,'[')
    endwhile
  br_position=strpos(calculate_par_str,']')
  while br_position ne -1 do begin
    calculate_par_str=strmid(calculate_par_str,0,br_position)+strmid(calculate_par_str,br_position+1)
    br_position=strpos(calculate_par_str,']')
    endwhile
;  print,calculate_par_str
  t=execute("calcP"+SSI(I)+"="+calculate_par_str)
  t=execute("if min(calcP"+SSI(I)+") - max(calcP"+ssi(i)+") eq 0 then skipcalcp=1")
  
;  if titles[n_elements(guess)+i] eq 'Broad to total flux fraction' then begin
;    index=where(calcp1 lt -3)
;    forprint,calcp1[index],ch12[index],ch10[index],ch2[index],ch0[index]
;    endif
;  if titles[n_elements(guess)+i] eq 'Broad to total flux fraction' then extra_str=',mininput=-1.0,maxinput=1.0' else $
;  if titles[n_elements(guess)+i] eq '6716/6731 broad' then extra_str=',mininput=-0.5,maxinput=2.8' else 
;  
;  if titles[n_elements(guess)+i] eq 'Broad to total flux fraction' then print,'MESSING WITH THINGS'
  extra_str=''
  
  if skipcalcp eq 1 then goto,skipcalcp_location
  if keyword_set(makeplots) THEN t=execute("cghistoplot,calcP"+SSI(I)+$
    ",nbins="+ssi(n_bins)+",title=titles[n_elements(guess)+i],histdata=histdata,locations=locations"+extra_str) $
      ELSE $
        t=execute("histdata=histogram(calcP"+SSI(I)+",nbins="+ssi(n_bins)+",locations=locations)",1,1)
  
  
  cup_out=mpmcmc_cup(locations,histdata)
  for ii=0,6 do cup_pars[i+n_elements(guess),ii]=cup_out[ii]
  if keyword_set(makeplots) then for ii=0,6 do cgplot,[cup_out[ii],cup_out[ii]],minmax(histdata),/overplot
  
  
  ;calculate initial guess and initial conditions.
  bestval=calc_mom_mpmcmcfun(locations,histdata,1)
  sigval=calc_mom_mpmcmcfun(locations,histdata,2)
  pi2 =[{fixed:0, limited:[1,1], limits:[0.D,max(histdata)*1.5]},$ ;peak value
    {fixed:0, limited:[0,0], limits:[0.D,0.d]},$ ;peak centroid
    {fixed:0, limited:[0,0], limits:[0.D,0.D]}] ;sigma
  estimates=[max(histdata),bestval,sigval]

  ;do the fit
  dummy=MPFITPEAK(locations,histdata,$
    coeff,nterms=3,/gaussian,$
    estimates=estimates,parinfo=pi2,status=status,chisq=chisq)
  if status lt 1 then begin
    ;plot the data.
    IF KEYWORD_SET(makeplots) THEN BEGIN
      cgplot,locations,gaussian(locations,estimates),/overplot,linestyle=2
      cgplot,locations,dummy,/overplot
      cgtext,0.2,0.9,ssf(bestval),/normal
      cgtext,0.2,0.85,ssf(sigval),/normal
      ENDIF
    return_value=[return_value,bestval]
    return_value_errs=[return_value_errs,sigval]
    gauss_fit_worked=[gauss_fit_worked,0]
    endif else begin;mpfitpeak status 1
      ;plot the data.
      IF KEYWORD_SET(makeplots) THEN BEGIN
        cgplot,locations,gaussian(locations,estimates),/overplot,linestyle=2
        cgplot,locations,dummy,/overplot
        cgtext,0.2,0.9,ssf(coeff[1]),/normal
        cgtext,0.2,0.85,ssf(coeff[2]),/normal
        cgtext,0.2,0.75,ssf(bestval),/normal
        cgtext,0.2,0.70,ssf(sigval),/normal
        ENDIF
      return_value=[return_value,coeff[1]]
      return_value_errs=[return_value_errs,coeff[2]]
      gauss_fit_worked=[gauss_fit_worked,1]
      endelse
  skipcalcp_location:
  if skipcalcp eq 1 then begin
    t=execute("return_value=[return_value,min(calcP"+SSI(I)+")]")
    return_value_errs=[return_value_errs,0.0]
    gauss_fit_worked=[gauss_fit_worked,0]
    endif ;skipcalcp
    

  if keyword_set(plotchain) then begin
    if n_iter gt 1000002 then $
      t=execute("cgplot,calcP"+SSI(I)+"[0:100000],/ynozero") $
        else t=execute("cgplot,calcP"+SSI(I)+",/ynozero")
      t=execute("cgplot,[0,1e7],[cup_out[0],cup_out[0]],/overplot")
      t=execute("cgplot,[0,1e7],[cup_out[1],cup_out[1]],/overplot")
      t=execute("cgplot,[0,1e7],[cup_out[2],cup_out[2]],/overplot")
      endif
  
      ;;CONTOURS.
      ;plot all contours.
      if keyword_set(makecontours) then begin
        for j=0,n_elements(guess)-1 do begin
          skipcontour=0
          t=execute("if min(CH"+SSI(J)+") - max(CH"+ssi(J)+") eq 0 then skipcontour=1")
          
          if parinfo[j].tied ne '' and parinfo[j].tied ne ' ' then skipcontour=1
          if parinfo[j].fixed eq 1 then skipcontour=1
          
          
          if skipcontour eq 1 then goto,skipcontour_location2
          t=execute("xr=minmax(calcP"+SSI(I)+")")
          t=execute("yr=minmax(CH"+SSI(j)+")")
;            ENDIF
          t=execute("binx=(max(calcP"+SSI(I)+")-min(calcP"+SSI(I)+"))/15.0")
          t=execute("biny=(max(ch"+SSI(j)+")-min(ch"+SSI(j)+"))/15.0")
          t=execute("density=hist_2d(calcP"+SSI(I)+",ch"+SSI(j)+",bin1=binx,bin2=biny,min1=min(calcP"+SSI(I)+"),max1=max(calcP"+SSI(I)+$
            "),min2=min(ch"+SSI(j)+"),max2=max(ch"+SSI(j)+"))")
          
          ;calcluate bin locations
          t=execute("x_locations=[min (calcP"+SSI(I)+")]")
          while n_elements(x_locations) lt n_elements(density[*,0]) do x_locations=[x_locations,x_locations[-1]+binx]
          t=execute("y_locations=[min (ch"+SSI(J)+")]")
          while n_elements(y_locations) lt n_elements(density[0,*]) do y_locations=[y_locations,y_locations[-1]+biny]
          x_locations=x_locations+binx/2
          y_locations=y_locations+biny/2
          contour,density,x_locations,y_locations,/ys,/xs,xr=xr,yr=yr,xtitle=titles[i+n_elements(guess)],ytitle=titles[j]
;          t=execute("cgplot,CH"+SSI(I)+",CH"+SSI(j)+",xtitle=titles[i],ytitle=titles[j]")
            skipcontour_location2:
            
          endfor;contour j
        endif;contour keyword.  
  
  
  
  
  
  
  endfor  ;calculated_params
  


  
  
  
  
  
  
  
  
if ~keyword_set(silent) then print,'time taken:',abs(TIMEvar-systime(/seconds)),' sec'
PERROR=return_value_errs
return,return_value
end
