function haganvol=Hagan_vol(alpha,nu,beta,rho,strike,T)
forward=1.0;
if strike==forward
    part_1 = (1.0 - beta)^2.0*alpha^2.0/(24.0*forward^(2.0 - 2.0*beta));
    part_2 = rho*beta*alpha*nu/(4.0*forward^(1.0 - beta));
    part_3 = (2.0 - 3.0*rho^2)*nu^2.0/24.0;
    haganvol=(alpha/forward^(1 - beta))*(1 + (part_1 + part_2 + part_3)*T );
else
    logfK = log(forward/strike);
    fkbpow = (forward*strike)^((1.0 - beta)/2.0);
    z = nu*fkbpow*logfK/alpha;
    xz = log((sqrt(1.0 - 2.0*rho*z + z^2.0 ) + z - rho)/(1.0-rho)) ;       
    part_1 = ((1.0-beta)^2.0)*(alpha^2.0)/(24.0*fkbpow^2.0);
    part_2 = (rho*beta*nu*alpha)/(4.0*fkbpow);
    part_3 = (2.0-3.0*rho^2)*nu^2.0/24.0;
    part_4 = ((1.0-beta)^2)*(logfK^2)/24.0;
    part_5 = ((1.0-beta)^4)*(logfK^4)/1920.0;
    haganvol = (alpha*z*(1 + (part_1 + part_2 + part_3)*T ))/(fkbpow*xz*(1 + part_4 + part_5 ));
end
end
    
    
    
    
