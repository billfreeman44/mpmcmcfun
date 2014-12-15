;returns log likelihood.
function likelih_mpmcmcfun,x_in,y_in,y_err,y_out
;return,total((-1.0)*(y_out-y_in)^2)>(-250)
denom=2.d*y_err*y_err
diff=y_in-y_out
x= total((-( diff*diff )/denom))
return,x>(-550)
end