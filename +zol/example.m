function [pts,val,info] = example(CAS)

% E = -1 (left) and F = +1 (right)
bnd         = [-1 1]; 
info.z4     = [];
info.z4x    = [];
switch CAS
    case '1a'
        Xlim    = 2.5*[-1 1];
        Ylim    = 1.5*[-1 1];
        theta   = linspace(0,2*pi,201);
        S       = exp(1i*theta);
        E       = -1+.5*S;
        F       = 1+.5*S;
        % Optimal Zolotarev
        v = ver;
        if any(strcmp('Symbolic Math Toolbox', {v.Name}))
            for i = 1:30
                syms x
                z           = zol.ZolOpt_1a(1,1/2,i);
                z4          = z(x);
                [num,den]   = numden(z4);
                num         = sym2poly(num); 
                den         = sym2poly(den);
                den_norm    = den(1);
                num         = num/den_norm; 
                den         = den/den_norm;
                info.z4{i}  = {[0;num(:)] den(:)};
                info.z4x{i} = z;
            end
        end
    case '1a2'
        Xlim    = 2.5*[-1 1];
        Ylim    = 1.5*[-1 1];
        S       = exp(2i*pi*(1:80)'/80);
        E       = -1.25+S;
        F       = 1.25+S;
        % Optimal Zolotarev
        v = ver;
        if any(strcmp('Symbolic Math Toolbox', {v.Name}))
            for i = 1:20
                syms x
                z           = zol.ZolOpt_1a(1,1/2,i);
                z4          = z(x);
                [num,den]   = numden(z4);
                num         = sym2poly(num); 
                den         = sym2poly(den);
                den_norm    = den(1);
                num         = num/den_norm; 
                den         = den/den_norm;
                info.z4{i}  = {[0;num(:)] den(:)};
                info.z4x{i} = z;
            end
        end
    case '1b'
        Xlim    = 2.5*[-1 1];
        Ylim    = 1.5*[-1 1];
        E       = zol.chebspace(-1.5,-.5,200);
        F       = zol.chebspace(.5,1.5,200);
        % Optimal Zolotarev
        v = ver;
        if any(strcmp('Symbolic Math Toolbox', {v.Name}))
            for i = 2:2:20
                syms x
                z           = zol.ZolOpt_1b_improper(1/2,3/2,i);
                z4          = z(x);
                [num,den]   = numden(z4);
                num         = sym2poly(num); 
                den         = sym2poly(den);
                den_norm    = den(1);
                num         = num/den_norm;
                den         = den/den_norm;
                info.z4{i}  = {[num(:)] den(:)};
                info.z4x{i} = z;
            end
        end
    case '1c'
        Xlim    = 2.5*[-1 1];
        Ylim    = 1.5*[-1 1];
        theta   = linspace(0,2*pi-1e-3,200);
        S       = exp(1i*theta);
        E       = -1+1i*zol.chebspace(-.75,.75,200);
        F       =  1+(.2*real(S)+1i*imag(S))/sqrt(1i);
    case '1d'
        Xlim    = 2.5*[-1 1];
        Ylim    = 1.5*[-1 1];
        theta   = linspace(pi/2,-pi/2,100);
        T       = .5*exp(1i*theta);
        TT      = exp(1i*theta);
        T1      = [T-.5+.5i -T-.5-.5i -TT-.5];
        E       = T1;
        F       = -T1;
    case '1e'
        Xlim    = 2.5*[-1 1];
        Ylim    = 1.5*[-1 1];
        theta   = linspace(-pi/2,pi/2,101);
        T       = exp(1i*theta);
        F       = 1-.74*T;
        Et      = zol.chebspace(-1.5,-.5,100)+.5*1i;
        Eb      = zol.chebspace(-1.5,-.5,100)-.5*1i;
        El      = -1.5+1i*zol.chebspace(-.5,.5,100);
        E       = [Et El Eb];
    case '1f'
        Xlim    = 2.5*[-1 1];
        Ylim    = 1.5*[-1 1];
        E       = 1-logspace(0,5,200);
        F       = zol.chebspace(1,2,200);
    case '2a'
        Xlim    = 2.5*[-1 1];
        Ylim    = 2*[-1 1];
        theta   = linspace(0,2*pi,200);
        S       = exp(1i*theta);
        E       = -1+.5*S;
        F1      = .8+.3*S+1i*.6;
        F2      = .8+.3*S-1i*.6;
        F       = [F1 F2];
    case '2b'
        Xlim    = 2.5*[-1 1];
        Ylim    = 2*[-1 1];
        F       = zol.chebspace(-.5,.5,100);
        E1      = zol.chebspace(-2,-1,100);
        E2      = zol.chebspace(1,2,100);
        E       = [E1 E2];
    case '2c'
        Xlim    = 2.5*[-1 1];
        Ylim    = 2*[-1 1];
        X1      =    zol.chebspace(-.5,.5,100);
        X2      = 1i*zol.chebspace(-.5,.5,100);
        X       = [X1 X2];
        E1      = X+1+1i;
        E2      = X-1-1i;
        E       = [E1 E2];
        F1      = X-1+1i;
        F2      = X*exp(1i*pi/4)+1-1i;
        F       = [F1 F2];
    case '2d'
        Xlim    = 2.5*[-1 1];
        Ylim    = 2*[-1 1];
        T1      = zol.chebspace(0,.5,100);
        T2      = zol.chebspace(0,.25,100)+1i*zol.chebspace(0,sqrt(.5^2-.25^2),100);
        T3      = zol.chebspace(.25,.5,100)-1i*(zol.chebspace(0,sqrt(.5^2-.25^2),100) - sqrt(.5^2-.25^2));
        T       = [T1 T2 T3]-.25-1i*.25;
        Tr      = T*exp(1i*pi/6);
        E       = [Tr-1 Tr-1+1i-1/sqrt(3) Tr-1-1i+1/sqrt(3)];
        F       = [T+1  T+1+1i            T+1-1i];
    case '3a'
        N       = 200; % W/T paper 
        %N       = 1000; 
        Xlim    = 2*[-1 1];
        Ylim    = 1.5*[-1 1];
        theta   = linspace(0,2*pi,N);
        %theta   = linspace(0,2*pi,2000);
        S       = exp(1i*theta);
        E       = .2+.5*S;
        F       = S;
    case '3b'
        Xlim    = 2*[-1 1];
        Ylim    = 1.5*[-1 1];
        theta   = linspace(0,2*pi,200);
        S       = exp(1i*theta);
        C1      = zol.chebspace(-.5,.5,100)-1i*.5;
        C2      = zol.chebspace(-.5,.5,100)+1i*.5;
        C3      = 1i*zol.chebspace(-.5,.5,100)-.5;
        C4      = 1i*zol.chebspace(-.5,.5,100)+.5;
        C       = [C1 C2 C3 C4];
        E       = C*exp(1i*pi/8);
        F       = 1.5*real(S)+1i*imag(S);
    case '3c'
        Xlim    = 2*[-1 1];
        Ylim    = 1.5*[-1 1];
        theta   = linspace(0,2*pi,200);
        S       = exp(1i*theta);
        E       = exp(.2i)*(.4*real(S)+.7i*imag(S));
        F       = S;
    case '3d'
        Xlim    = 2*[-1 1];
        Ylim    = 1.5*[-1 1];
        theta   = linspace(0,2*pi,200);
        S       = exp(1i*theta);
        E       = [.1+.5i+.3*S .1-.5i+.3*S];
        F       = S;
    case '7'
        Xlim    = 1.5*[-1 1];
        Ylim    = 1.5*[-1 1];
        C1      = 1/4+1i*zol.chebspace(-1,1,100);
        C2      = zol.chebspace(1/4,1,50)+1i;
        C3      = 1+1i*zol.chebspace(-1,1,100);
        C4      = zol.chebspace(1/4,1,50)-1i;
        C       = [C1 C2 C3 C4];
        E       = -C;
        F       = C;
    case 'spiral1'
        Xlim    = 55*[-1 1];
        Ylim    = 55*[-1 1];
        t       = 1:200; 
        u       = .0265; 
        r0      = .1;
        r       = r0 + u*10*t;
        omega   = .005;
        phi0    = 3*pi/2;
        phi     = -omega*10*t+phi0;
        x       = r .* cos(phi);
        y       = r .* sin(phi);
        E       = x + 1i*y;
        F       = -x - 1i*y;
    case 'pm2'
        N       = 120;
        Xlim    = 2.5*[-1 1];
        Ylim    = 1.5*[-1 1];
        theta   = linspace(pi/4,2*pi-pi/5,2*N);
        S       = exp(1i*theta);
        E       = -1+S;
        B1      = linspace(-1,real(E(1)),N/2+1)  +1i*linspace(0,imag(E(1)),N/2+1);
        B2      = linspace(-1,real(E(end)),N/2+2)+1i*linspace(0,imag(E(end)),N/2+2);
        theta   = linspace(0,2*pi-1e-5,40);
        O       = -1+1i*.75+.1*exp(1i*theta);
        E       = [E B1(1:end-1) B2(2:end-1) O];
        %
        theta   = linspace(0,2*pi,numel(E));
        S       = exp(1i*theta);
        F       = 1+S;
    case 'spiral2'
        Xlim    = 55*[-1 1];
        Ylim    = 55*[-1 1];
        t = 1:200; u = .0265; r0 = .1;
        r = r0 + u*10*t;
        omega = .005;
        phi0 = 3*pi/2;
        phi = -omega*10*t+phi0;
        x = r .* cos(1.2*phi);
        y = r .* sin(1.2*phi);
        E = x + 1i*y;
        F = -x - 1i*y;% n = 2000;
end

%%% Data
la          = E(:); 
n1          = length(la);
mu          = F(:); 
n2          = length(mu);
W(1:n1,1)   = min(bnd);
V(1:n2,1)   = max(bnd);
pts         = [la; mu];
val         = [W; V];

% 
info.E      = E;
info.F      = F;
info.bnd    = bnd; % [E F]
info.Xlim   = Xlim;
info.Ylim   = Ylim;


% %%
% function x = zol.chebspace(a,b,n)
% 
%     x = zeros(1,n);
%     for k = 0:n-1
%         x(1,k+1) = cos((k+.5)*pi/n);
%     end
% 
%     x = .5*(a+b)+.5*(b-a)*x;
% 
% end