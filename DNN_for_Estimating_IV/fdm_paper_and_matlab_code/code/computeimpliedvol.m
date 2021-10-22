function [MCvol,FDMvol]=computeimpliedvol(strike,tenor,value,fdm_value)
forward=1.0; T=tenor;
param_fdm=[forward,strike,T,fdm_value]; param_mc=[forward,strike,T,value];
fun1=@(vol)optim(vol,param_fdm); fun2=@(vol)optim(vol,param_mc);
init_guess=0.2;
MCvol=fminsearch(fun2,init_guess);FDMvol=fminsearch(fun1,init_guess);
end


function error=optim(vol,param)
forward=param(1);strike=param(2);T=param(3);price=param(4);
error=(price-BScall(forward,strike,T,vol))^2;
end

function call=BScall(forward,strike,T,sigma)
d1=(log(forward/strike) + sigma^2/2*T)/(sigma * sqrt(T));
d2= d1- sigma* sqrt(T);
call=forward * normcdf(d1,0,1) - strike * normcdf(d2,0,1);
end
