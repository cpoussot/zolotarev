function [Z4,z4,sigma] = ZolOpt_1b_improper(a,b,r)
    d           = r/2; 
    K           = ellipke(1-(a/b)^2);
    [sn,cn,dn]  = ellipj((0:2*d+1)*K/(2*d+1),1-(a/b)^2);
    c           = a^2*sn.^2./cn.^2;
    zero        = -c(3:2:end-1);
    pole        = -c(2:2:end-2);
    extr        = a^(2)*dn.^(-2);
    
    % Determine scaling factor D.
    rr  = 1;
    for j = 1:d
        rr = rr.*(extr-zero(j))./(extr-pole(j));
    end
    rr  = rr.*sqrt(extr);
    D   = 2/(min(rr)+max(rr));
    
    % now D*r is relative approximation to 1/sqrt(x) on [0,b^2]
    % now use sign(x)=x/sqrt(x^2)
    zer(1:2:2*d)    = sqrt(zero);
    zer(2:2:2*d)    = -sqrt(zero);
    pol(1:2:2*d)    = sqrt(pole);
    pol(2:2:2*d)    = -sqrt(pole);
    zer(2*d+1)      = 0;
    
    % Zolotarev optimal approximant as function handle
    Z4      = @(s) feval(@vhandle,s,pol,zer,D);
    z4      = z4_nd(pol,zer,D);
    sigma   = D;

    % % Zolotarev optimal approximant as function handle
    % z               = @(s) feval(@vhandle,s,pol,zer,D); 
    % [sn2,cn2,dn2]   = ellipj((0:4*d+2)*K/(4*d+2),1-(a/b)^2);
    % c2              = sqrt(a^2+(b^2-a^2)*cn2.^2);
    % intp            = fliplr(c2(2:2:end));
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

function z4 = z4_nd(pol,zer,D)
    num = 1;
    for jj = 1:length(zer)
        num = conv(num,[1 -zer(jj)]);
    end
    den = 1;
    for jj = 1:length(pol)
        den = conv(den,[1 -pol(jj)]);
    end
    z4 = tf(D*num,den);
end