;returns log likelihood.
function likelih_mpmcmcfun,x_in,y_in,y_err,y_out
;http://www.r-bloggers.com/maximum-likelihood/
x= total((-( (y_in-y_out)*(y_in-y_out) )/(2.0*y_err*y_err)))
return,x>(-687) ;405?
end