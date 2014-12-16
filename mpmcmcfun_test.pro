pro mpmcmcfun_test

;create fake data
x=findgen(100)/10.0
nel=n_elements(x)
errs=replicate(0.1d,nel)
y=replicate(0.0d,nel)
x=findgen(100)/10.0
n_iter=30000

y=gaussian(x,[1.1,5.2,1.6])
y=y +randomn(seed,nel)*errs

 p=[1.0,5.0,1.5]
   pi =[{fixed:0, limited:[1,1], limits:[0.D,1.2d], tied:''},$ ;peak value  
        {fixed:0, limited:[1,1], limits:[4.0d,5.5d], tied:''},$ ;peak centroid 1
        {fixed:0, limited:[1,1], limits:[1.00D,2.00D], tied:''}] ;sigma 2]; LINEAR TERM 3
        
  result1=mpfitfun('gaussian',x,y,errs,p,parinfo=pi,/QUIET,status=status,PERROR=PERROR1)
;  stop
  ps_start,'~/mcmc_test.ps'
  result2=mpmcmcfun('gaussian',x,y,errs,result1,perror1,n_iter,parinfo=pi,status=status,PERROR=PERROR2,/makeplots,/makecontours,$
    gauss_fit_worked=gauss_fit_worked)
  cgplot,x,y,title='best fit using mpfitfun'
  cgplot,x,gaussian(x,result1),/overplot,color='red'
  cgplot,x,y,title='best fit using mpmcmcfun'
  cgplot,x,gaussian(x,result2),/overplot,color='red'
  ps_end

print
print,result1
print,perror1
print
print,result2
print,perror2


end