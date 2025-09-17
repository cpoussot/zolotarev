function [z,pol,zer,sigma] = ZolOpt_1a(alpha,rho,r)

%   7 July 2025
%   Author: I.V. Gosea


zer = ones(1,r)*sqrt(alpha^2-rho^2);
pol = -ones(1,r)*sqrt(alpha^2-rho^2);

sigma = ((1-zer(1))/(1+zer(1)))^r;

z = @(s) feval(@vhandle,s,pol,zer,sigma); 

end

function z = vhandle(s,pol,zer,sigma)

Num = 1;

for jj = 1:length(zer)
    Num = Num*(s-zer(jj));
end

Den = 1;

for jj = 1:length(pol)
    Den = Den*(s-pol(jj));
end

% for problem Z3
zinit = sigma^(1/2)*Num/Den;

% for problem Z4
z = -((1-sigma)/(1+sigma))*((zinit-sqrt(sigma))/(zinit+sqrt(sigma)));

end

