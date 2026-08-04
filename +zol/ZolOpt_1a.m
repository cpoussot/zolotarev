function [Z4,z4,sigma] = ZolOpt_1a(alpha,rho,r)
    zer     = ones(1,r)*sqrt(alpha^2-rho^2);
    pol     = -ones(1,r)*sqrt(alpha^2-rho^2);
    sigma   = ((1-zer(1))/(1+zer(1)))^r;
    Z4      = @(s) feval(@vhandle,s,pol,zer,sigma);
    z4      = z4_nd(pol,zer,sigma);
end

function Z4 = vhandle(s,pol,zer,sigma)
    Num = 1;
    for jj = 1:length(zer)
        Num = Num*(s-zer(jj));
    end
    Den = 1;
    for jj = 1:length(pol)
        Den = Den*(s-pol(jj));
    end
    % for problem Z3
    Z3  = sigma^(1/2)*Num/Den;
    % for problem Z4
    Z4  = -((1-sigma)/(1+sigma))*((Z3-sqrt(sigma))/(Z3+sqrt(sigma)));
end

function z4 = z4_nd(pol,zer,sigma)
    num = 1;
    for jj = 1:length(zer)
        num = conv(num,[1 -zer(jj)]);
    end
    den = 1;
    for jj = 1:length(pol)
        den = conv(den,[1 -pol(jj)]);
    end
    % for problem Z3
    z3_num = sigma^(1/2)*num;
    z3_den = den;
    % for problem Z4
    scaling = -((1-sigma)/(1+sigma));
    z4_N    = z3_num-sqrt(sigma)*z3_den;
    z4_D    = z3_num+sqrt(sigma)*z3_den;
    z4      = tf(scaling*z4_N,z4_D);
end

