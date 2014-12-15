;1: x location of y max
;2: approximate 1 sigma location based on literal FWHM
;   if FWHM cannot be determined then it returns
;   (max(x)-min(x))*0.25
;3: supposedly the first moment
;4: supposedly the second moment.
function calc_mom_mpmcmcfun,xfit,yfit,num
    if num eq 1 then begin 
      t=max(yfit,ind)
;      if abs(xfit[ind]) lt abs(xfit[0]-xfit[1]) then return,xfit[ind]+abs(xfit[0]-xfit[1])
      return,xfit[ind]
      endif
    if num eq 2 then begin
      t=max(yfit,ind)
      index=where(yfit lt max(yfit)*0.5,ct)
      if ct lt 2 then return,abs((max(xfit)-min(xfit))*0.25)
;      if ct lt 2 then message,'abs((max(xfit)-min(xfit))*0.25)'
      t=max(yfit[index],ind2)
      return,abs(xfit[ind] - xfit[index[ind2]])*2.0/2.35
      endif
    if num eq 3 then return,total(xfit*yfit)/total(yfit)
    if num eq 4 then return,sqrt(abs(total(xfit*xfit*yfit)/total(yfit) - $
      (total(xfit*yfit)/total(yfit))*(total(xfit*yfit)/total(yfit))));2.355*
message,'calc mom mcmc error, num needed and should be 1 or 2. (or 3 or 4 i guess)'
return,0.0
end
