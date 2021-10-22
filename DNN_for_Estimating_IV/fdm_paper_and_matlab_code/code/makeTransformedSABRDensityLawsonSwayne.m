function [P, PL, PR, zm, zmin, zmax, h] = makeTransformedSABRDensityLawsonSwayne(alpha, beta, nu, rho, forward, T, N, timesteps, nd)
[zmin, zmax] = computeBoundaries(alpha, beta, nu, rho, forward, T, nd); 
J = N-2; h0 = (zmax-zmin)/(J); j0 = int32((0-zmin)/h0);
h = (0-zmin)/(double(j0)-0.5);
z = (0:(J+1))*h + zmin; zmax = z(J+2); zm = z - 0.5*h;
ym = Y(alpha,nu,rho,zm);
ymax= Y(alpha,nu,rho,zmax); ymin=Y(alpha,nu,rho,zmin);
Fm=F(forward, beta, ym);
Fmax=F(forward, beta, ymax); Fmin = F(forward, beta, ymin);
Fm(1) = 2*Fmin-Fm(2); Fm(J+2)= 2*Fmax -Fm(J+1);
Cm = C(alpha, beta, rho, nu, ym, Fm); Cm(1) = Cm(2); Cm(J+2) = Cm(J+1);
Gammam = G(forward, beta, Fm, j0); 
dt = T/timesteps;
b = 1- sqrt(2)/2; % lawson Swayne param
dt1 = dt * b; dt2 = dt * (1-2*b); Em = ones(1,J+2);
Emdt1= exp(rho*nu*alpha*Gammam*dt1); Emdt1(1)=Emdt1(2); Emdt1(J+2)=Emdt1(J+1);
Emdt2= exp(rho*nu*alpha*Gammam*dt2); Emdt2(1)=Emdt2(2); Emdt2(J+2)=Emdt2(J+1);
PL = 0.0; PR = 0.0; P=zeros(J+2,1); P(j0+1,1)=1.0/h;

for t = 1:timesteps
    Em=Em.*Emdt1;[P1,PL1,PR1]=solveStep(Fm, Cm, Em, dt1,h,P,PL,PR);
    Em=Em.*Emdt1;[P2,PL2,PR2]=solveStep(Fm, Cm, Em, dt1,h,P1,PL1,PR1);% Emdt1 or EMdt2
    P=(sqrt(2)+1)*P2 -sqrt(2)*P1;
    PL=(sqrt(2)+1)*PL2 -sqrt(2)*PL1;
    PR=(sqrt(2)+1)*PR2 -sqrt(2)*PR1;
    Em=Em.* Emdt2;
end
end

function [P, PL, PR] = solveStep(Fm, Cm, Em, dt, h, P, PL, PR)
frac = dt/(2*h); M = length(P);
B(2:M-1) = 1.0 + frac*(Cm(2:M-1).*Em(2:M-1).*(1./(Fm(3:M)-Fm(2:M-1))+1./(Fm(2:M-1)-Fm(1:M-2))));
C(2:M-1) = -frac* Cm(3:M).*Em(3:M)./(Fm(3:M)-Fm(2:M-1));
A(1:M-2) = -frac* Cm(1:M-2).*Em(1:M-2)./(Fm(2:M-1)-Fm(1:M-2));
B(1) = Cm(1)/(Fm(2)-Fm(1))*Em(1); 
C(1) = Cm(2)/(Fm(2)-Fm(1))*Em(2);
B(M) = Cm(M)/(Fm(M)-Fm(M-1))*Em(M); 
A(M-1) = Cm(M-1)/(Fm(M)-Fm(M-1))*Em(M-1); 
tri = diag(sparse(B))+diag(sparse(A),-1)+diag(sparse(C),1);
%tri(2,1:end);
P(1) = 0; P(M) = 0;
P = tri\P;
PL = PL + dt*Cm(2)/(Fm(2)-Fm(1))*Em(2)*P(2);
PR = PR + dt*Cm(M-1)/(Fm(M)-Fm(M-1))*Em(M-1)*P(M-1);
end

function[zmin, zmax] = computeBoundaries(alpha, beta, nu, rho, forward, T, nd)
zmin = -nd*sqrt(T); zmax = -zmin; 
end

function F = F(forward , beta , ym)
u = sign(forward)*abs(forward)^(1-beta)+(1-beta)*ym; F = sign(u).*abs(u).^(1/(1-beta));
end 

function Cm = C(alpha, beta, rho, nu, ym, Fm)
Cm = sqrt(alpha^2+2*rho*alpha*nu*ym+nu^2*ym.^2).*abs(Fm).^(beta); 
end

function G = G(forward, beta, Fm, j0)
G = (abs(Fm).^beta-abs(forward)^beta)./(Fm-forward); 
G(j0+1) = sign(forward)*beta/abs(forward).^(1-beta);
end

% function y = yOfStrike(strike , forward , beta)
% y = (sign(strike)*abs(strike)^(1-beta)-sign(forward)*abs(forward)^(1-beta))/(1- beta);
% end
function Y = Y(alpha, nu, rho, zm)
Y = alpha/nu*(sinh(nu*zm)+rho*(cosh(nu*zm)-1));
end

