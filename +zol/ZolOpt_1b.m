function [z,pol,zer,D,intp] = ZolOpt_1b(a,b,r)

%   27 March 2019
%   Author: I.V. Gosea

d = r/2; 

K = ellipke(1-(a/b)^2);
%K = util_ellipkkp(1-(a/b)^2);
 
[sn,cn,dn] = ellipj((0:2*d)*K/(2*d),1-(a/b)^2);
%[sn,cn,dn] = util_ellipjc((0:2*r)*K/(2*r),1-(a/b)^2);

c    = a^2*sn.^2./cn.^2;
zero = -c(3:2:end-1);
pole = -c(2:2:end);

extr = a^(2)*dn.^(-2);

% Determine scaling factor D.
rr=1;
for j=1:d
    if j<d
        rr=rr.*(extr-zero(j))./(extr-pole(j));
    else
        rr=rr./(extr-pole(j));
    end
end
rr=rr.*sqrt(extr);
D=2/(min(rr)+max(rr));

% now D*r is relative approximation to 1/sqrt(x) on [0,b^2]
% now use sign(x)=x/sqrt(x^2)
zer(1:2:2*(d-1))=sqrt(zero);
zer(2:2:2*(d-1))=-sqrt(zero);
pol(1:2:2*d)=sqrt(pole);
pol(2:2:2*d)=-sqrt(pole);
zer(2*d-1)=0;

% Zolotarev optimal approximant as function handle
z = @(s) feval(@vhandle,s,pol,zer,D); 

[sn2,cn2,dn2] = ellipj((0:4*d)*K/(4*d),1-(a/b)^2);

%[sn2,cn2,dn2] = util_ellipjc((0:4*d)*K/(4*d),1-(a/b)^2);

 c2 =sqrt(a^2+(b^2-a^2)*cn2.^2);
 
 intp = fliplr(c2(2:2:end));

end

function z = vhandle(s,pol,zer,D)

Num = 1;

for jj = 1:length(zer)
    Num = Num*(s-zer(jj));
end

Den = 1;

for jj = 1:length(pol)
    Den = Den*(s-pol(jj));
end

z = D*Num/Den;

end

