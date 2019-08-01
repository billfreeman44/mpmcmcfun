function mpmcmc_cup,x,y

cup_out=replicate(-1d,7)

;best fit is peak value
  peakval=max(y,max_subscript)
  cup_out[0]=x[max_subscript]
  tot=total(y)
  
  ;calculate 1,2,3 sigma (68.27, 95.45, 99.73)
  high_ind=max_subscript
  low_ind=max_subscript
  current_total=peakval
  
  low_limit_reached=0
  high_limit_reached=0
  while current_total/tot lt 0.6827 do begin
    if low_ind-1 lt 0 then low_limit_reached=1
    if high_ind+1 ge n_elements(y) then high_limit_reached=1
    if high_limit_reached and low_limit_reached then goto,loopbail1
    
    
    if low_limit_reached then begin
      high_ind=high_ind+1
      current_total=current_total+y[high_ind]
      endif else $
        if high_limit_reached then begin
          low_ind=low_ind-1
          current_total=current_total+y[low_ind]
          endif else $
            if y[low_ind-1] gt y[high_ind+1] then begin
              low_ind=low_ind-1
              current_total=current_total+y[low_ind]
              endif else $
                if y[high_ind+1] ge y[low_ind-1] then begin
                  high_ind=high_ind+1
                  current_total=current_total+y[high_ind]
                  endif ;else goto,loopbail1

    endwhile ; 1 sigma
  loopbail1:
  cup_out[1]=x[high_ind]
  cup_out[2]=x[low_ind]
  
  
  low_limit_reached=0
  high_limit_reached=0
  while current_total/tot lt 0.9545 do begin
    if low_ind-1 lt 0 then low_limit_reached=1
    if high_ind+1 ge n_elements(y) then high_limit_reached=1
    if high_limit_reached and low_limit_reached then goto,loopbail2
    
    
    if low_limit_reached then begin
      high_ind=high_ind+1
      current_total=current_total+y[high_ind]
      endif else $
        if high_limit_reached then begin
          low_ind=low_ind-1
          current_total=current_total+y[low_ind]
          endif else $
            if y[low_ind-1] gt y[high_ind+1] then begin
              low_ind=low_ind-1
              current_total=current_total+y[low_ind]
              endif else $
                if y[high_ind+1] ge y[low_ind-1] then begin
                  high_ind=high_ind+1
                  current_total=current_total+y[high_ind]
                  endif ;else goto,loopbail1
    endwhile ; 1 sigma
  loopbail2:
  cup_out[3]=x[high_ind]
  cup_out[4]=x[low_ind] 
  
  
  low_limit_reached=0
  high_limit_reached=0
  while current_total/tot lt 0.99 do begin
    if low_ind-1 lt 0 then low_limit_reached=1
    if high_ind+1 ge n_elements(y) then high_limit_reached=1
    if high_limit_reached and low_limit_reached then goto,loopbail3
    
    
    if low_limit_reached then begin
      high_ind=high_ind+1
      current_total=current_total+y[high_ind]
      endif else $
        if high_limit_reached then begin
          low_ind=low_ind-1
          current_total=current_total+y[low_ind]
          endif else $
            if y[low_ind-1] gt y[high_ind+1] then begin
              low_ind=low_ind-1
              current_total=current_total+y[low_ind]
              endif else $
                if y[high_ind+1] ge y[low_ind-1] then begin
                  high_ind=high_ind+1
                  current_total=current_total+y[high_ind]
                  endif ;else goto,loopbail1
    endwhile ; 1 sigma
  loopbail3:
  cup_out[5]=x[high_ind]
  cup_out[6]=x[low_ind]  
  



return,cup_out
end