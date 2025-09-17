function [pts,val,info] = example(CAS)

% E = -1 (left) and F = +1 (right)
bnd     = [-1 1]; 
%bnd     = [0 1]; 
info.z4 = [];
switch CAS
    case '1a'
        Xlim    = 2.5*[-1 1];
        Ylim    = 1.5*[-1 1];
        theta   = linspace(0,2*pi,201);
        S       = exp(1i*theta);
        E       = -1+.5*S;
        F       = 1+.5*S;
        % Optimal Zolotarev
        % info.z4{2}  = {[0 1.71429 0] [1 0 .75]};
        % info.z4{3}  = {[0 2.59615 0 .649038] [1 0 2.25 0]};
        % info.z4{4}  = {[0 3.46392 0 2.59794 0] [1 0 4.5 0 .5625]};
        % % 
        for i = 1:10
            [z,pol,zer,sigma] = zol.ZolOpt_1a(1,1/2,i);
            syms x
            %zz = zz(x);
            % DEN = 1;
            % for ii = 1:length(pol)
            %     DEN = DEN * (x-pol(ii));
            % end
            % NUM = 1; 
            % for ii = 1:length(zer)
            %     NUM = NUM * (x-zer(ii));
            % end
            % z3          = sigma^(1/2)*NUM/DEN;
            % z4          = -((1-sigma)/(1+sigma))*((z3-sqrt(sigma))/(z3+sqrt(sigma)));
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
    case '1b'
        Xlim    = 2.5*[-1 1];
        Ylim    = 1.5*[-1 1];
        E       = zol.chebspace(-1.5,-.5,200);
        F       = zol.chebspace(.5,1.5,200);
        % Optimal Zolotarev
        % info.z4{2}  = {[0 1.71429 0] [1 0 .75]};
        % info.z4{3}  = {[0 2.59615 0 .649038] [1 0 2.25 0]};
        % info.z4{4}  = {[0 3.46392 0 2.59794 0] [1 0 4.5 0 .5625]};
        % % 
        for i = 2:2:10
            [z,pol,zer,sigma,~] = zol.ZolOpt_1b_improper(1/2,3/2,i);
            %[z,pol,zer,sigma,~] = zol.ZolOpt_1b(1/5,3/2,i);
            syms x
            % DEN = 1;
            % for ii = 1:length(pol)
            %     DEN = DEN * (x-pol(ii));
            % end
            % NUM = 1; 
            % for ii = 1:length(zer)
            %     NUM = NUM * (x-zer(ii));
            % end
            % z3          = sigma^(1/2)*NUM/DEN;
            % z4          = -((1-sigma)/(1+sigma))*((z3-sqrt(sigma))/(z3+sqrt(sigma)));
            z4          = z(x);
            [num,den]   = numden(z4);
            num         = sym2poly(num); 
            den         = sym2poly(den);
            den_norm    = den(1);
            num         = num/den_norm;
            den         = den/den_norm;
            %info.z4{i}  = {[0 num] den};
            info.z4{i}  = {[num(:)] den(:)};
            info.z4x{i} = z;
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
    case 'pm'
        N       = 120;
        Xlim    = 2.5*[-1 1];
        Ylim    = 1.5*[-1 1];
        theta   = linspace(pi/4,2*pi-pi/5,2*N);
        S       = exp(1i*theta);
        E       = -1+S;
        B1      = linspace(-1,real(E(1)),N/2+1)  +1i*linspace(0,imag(E(1)),N/2+1);
        B2      = linspace(-1,real(E(end)),N/2+2)+1i*linspace(0,imag(E(end)),N/2+2);
        theta   = linspace(0,2*pi,37);
        O       = -1+1i*.75+.1*exp(1i*theta-1e-12);
        E       = [E B1(1:end-1) B2(2:end-1) O];
        nE      = length(E);
        %
        % theta   = linspace(pi/2,-pi/2,N);
        % T       = .5*exp(1i*theta);
        % TT      = exp(1i*theta);
        % T1      = [T-.5+.5i -T-.5-.5i -TT-.5];
        % F       = -.5*T1*exp(1i*pi/5);
        A1  = linspace(0,.4,floor(nE/9))+1i*linspace(0,1,floor(numel(E)/9));
        A2  = linspace(.41,.8,floor(numel(E)/9))+1i*linspace(1,0,floor(numel(E)/9));
        A3  = linspace(.16,.64,floor(numel(E)/9))+1i*.4;
        %A3  = linspace(.18,.62,floor(numel(E)/9))+1i*.4;
        A   = .8*([A1 A2 A3]-1i*.5);
        F   = .8*[A A+.8 A+1.6]+.5+.5i;
        %F   = .8*[A exp(-1i*pi/10)*A+.9 exp(-1i*2*pi/10)*A+1.8]+.5;
        %F   = [.8*A+1i .8*A+.5 .8*A-1i]+.5;
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
    case 'ts' % aca (triangles and squares)
        N       = 50;
        Xlim    = 3*[-1 1];
        Ylim    = 3*[-1 1];
        p1=1+1i*1; p2=1+1i*2; p3=2+1i*2; p4=2+1i;
        p12=linspace(p1,p2,25);p23=linspace(p2,p3,25);
        p34=linspace(p3,p4,25);p41=linspace(p4,p1,25);
        P=union(union(union(p12,p23),p34),p41);

        q1=1-1i; q2=2-1i; q3=3/2-1i*2;
        q12=linspace(q1,q2,N);q23=linspace(q2,q3,N);q31=linspace(q3,q1,N);
        Q=union(union(q12,q23),q31);
        
        r1=q1-3; r2=q2-3; r3=q3-3;
        r12=linspace(r1,r2,N);r23=linspace(r2,r3,N);r31=linspace(r3,r1,N);
        R=union(union(r12,r23),r31);
        
        m1=-1+1i; m2=-2+1i; m3=-3/2+1i*2;
        m12=linspace(m1,m2,N);m23=linspace(m2,m3,N);m31=linspace(m3,m1,N);
        M=union(union(m12,m23),m31);
        
        E   = union(P,R); 
        F   = union(Q,M);
    case 1
        N       = 500;
        Xlim    = 2.5*[-1 1];
        Ylim    = 1.5*[-1 1];
        theta   = linspace(0,2*pi-1e-3,100);
        S       = exp(1i*theta);
        F       = 1+exp(1i*pi/2)*(.2*real(S)+1i*imag(S))/sqrt(1i);
        % Polygon random
        pgon    = polyshape(rand(1,N),1*rand(1,N));
        E       = pgon.Vertices(:,1)+1i*pgon.Vertices(:,2); E(isnan(E)) = [];
        k       = boundary([real(E) imag(E)]);
        E       = E(k,:)*exp(1i*pi/2*rand(1))-1.1*max(real(E));
        E(1)    = []; % to ensure (la U mu) not repeating
    case 3  % aca Two line segments
        Xlim    = .3*[-1 1];
        Ylim    = .2*[-1 1];
        F       = 1i*([0.005:0.005:1]/10+0.01);
        F       = double(F); 
        syms s;
        dt = 0.005;
        z0 = subs(s+.1,s,[0:0.005:1-0.005]);
        E  = [1:200]/200+i*z0/10+.1;E=double(E);
        [n1,n2]=size(E);
        n=n2;
    case 4 % aca (two circles)
        Xlim    = 2.5*[-1 1];
        Ylim    = 1.5*[-1 1];
        dt=0.005; th=[0:dt:1-dt]*2*pi; [n1,n2]=size(th);n=n2;
        S=exp(1i*th);
        F=-.2+.2*1i+S/2;%F +1
        E=S;            %E -1 or 0
    case 5
        Xlim    = 2.5*[-1 1];
        Ylim    = 3.5*[-1 1];
        n0=100;
        E1=linspace(1,2,n0); E2=linspace(-3,3,n0)*1i;
        E3=E1+1i*3;E7=-E1+1i*3; 
        E4=E1-1i*3; E8=-E1-1i*3;
        E5=E2+1;E9=E2-1; E6=E2+2; E10=E2-2;
        E=union(union(union(E3,E4),E5),E6);
        F=union(union(union(E7,E8),E9),E10); 
        [n1,n2]=size(E);n=n2;
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