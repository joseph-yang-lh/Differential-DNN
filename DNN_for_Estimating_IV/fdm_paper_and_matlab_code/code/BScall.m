function call=BScall(forward,strike,T,sigma)
d1=(log(forward/strike) + sigma^2/2*T)/(sigma * sqrt(T));
d2= d1- sigma* sqrt(T);
call=forward * normcdf(d1,0,1) - strike * normcdf(d2,0,1);
end