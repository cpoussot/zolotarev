clearvars; close all; clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN: TO ADAPT BY USER
% Add MDSPACK
%%% MDSPACK
setenv('MDSHOME','/Users/charles/Documents/MDS')
addpath('/Users/charles/Documents/MDS/mdspack/MDSPACK/osx/v1.1.0/API/matlab/')
addpath('/Users/charles/Documents/MDS/mdspack/MDSPACK/osx/v1.1.0/bin')
% Add Zolotarev Loewner package
addpath('/Users/charles/Documents/GIT/zolotarev')
addpath('/Users/charles/Documents/GIT/lf')
% Add AAA package
addpath('/Users/charles/Documents/GIT/_others/chebfun')
%%% Choose case, order
CAS     = 'iir';
robj0   = 1e-12; % objective order (either integer > 1 or sigma threshold)
REALTF  = true;
SAVEIT  = false;
Fs      = 5;
Ts      = 1/Fs;
%%% END: TO ADAPT BY USER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Plot properties
set(groot,'DefaultFigurePosition', [200 150 1000 600]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',4)
set(groot,'defaultaxesfontsize',18)
set(groot,'defaultAxesTickLabelInterpreter','latex');  
list_facto  = fieldnames(get(groot,'factory'));index_interpreter = find(contains(list_facto,'Interpreter'));for i = 1:length(index_interpreter); set(groot, strrep(list_facto{index_interpreter(i)},'factory','default'),'latex'); end
mw          = 15;         % marker width
lev_contour = -30:1:0;    % contour plot
col         = parula(10); % 
col1        = col(5,:);   % Z4 poles
col2        = col(9,:);   % Z4 zeros

%%% Define Zolotarev topology
[pts,val,data]  = zol.example(CAS);

%%% Loewner approximation
% >> (Z4) rational approximation
tic
methodName      = '\textbf{LF}';
[la,mu,W,V]     = zol.example2data(pts,val,data,true,REALTF);
opt             = [];
opt.target      = robj0;
opt.real        = REALTF;
[h4,info]       = zol.loewner(la,mu,W,V,opt);

% >> (Z4) poles and zeros
h4poles         = eig(info.Ar,info.Er);
h4zeros         = eig([info.Ar info.Br;info.Cr 0],blkdiag(info.Er,0));
timeLOE         = toc;
robj            = info.r; 
% >> (Z3) rational approximation, from Z4->Z3
[h3,hp,hsig]    = zol.pb4_to_pb3(h4,pts,val);
% >> (Z3) poles and zeros
h3poles         = eig([info.Ar info.Br;-info.Cr (hp)],blkdiag(info.Er,0));
h3zeros         = eig([info.Ar info.Br; info.Cr (hp)],blkdiag(info.Er,0));
% >> (Z3) evaluate
xx              = linspace(data.Xlim(1),data.Xlim(2),101);
yy              = linspace(data.Ylim(1),data.Ylim(2),103);
[X,Y]           = meshgrid(xx,yy);
Zr3_loe         = zeros(numel(yy),numel(xx));
for i = 1:numel(xx)
    for j = 1:numel(yy)
        Zr3_loe(j,i) = h3(xx(i)+1i*yy(j));
    end
end

%%% Stable equivalent digital filter
Fz      = dss(info.Ar,info.Br,info.Cr,info.Dr,info.Er,Ts);
Fz      = zol.mirror(Fz);
h4      = @(z) Fz.c*((Fz.e*z-Fz.a)\Fz.b) + Fz.d;
eigF    = eig(Fz);

% figure
% nyquist(        Fz), hold on
% nyquist(    d2c(Fz,'tustin'))
% nyquist(c2d(d2c(Fz,'tustin'),1/Fs,'tustin'))

%
%%% pH (continuous time)
Fc      = d2c(Fz,'tustin');
h4c     = @(s) Fc.c*((Fc.e*s-Fc.a)\Fc.b) + Fc.d;
Ds      = norm(Fc,inf)*1.1;
D       = Fc.d;
Fcd     = Fc-D+Ds;
puls1   = logspace(-2,log10(pi*Fs),1e4);
[hloeph,info_loeph] = lf.passive2ph(Fcd);
hloeph = info_loeph.Hr;
hloeph.E = eye(length(hloeph.A));
hloeph.D = info_loeph.Hr.D+D-Ds;
hloeph = c2d(hloeph,Ts,'tustin');
hloeph = @(z) hloeph.c*((hloeph.e*z-hloeph.a)\hloeph.b) + hloeph.d;
%info_loeph.S = info_loeph.S +D-Ds;
% figure, hold on
% mdspack.nyquist(Fc,puls1)
% mdspack.nyquist(Fcd,puls1)
% mdspack.nyquist(hloeph,puls1,'--')
% %xlim([-1 1]), ylim([-1 1])

% %
% nip     = 100;
% puls    = logspace(-2,log10(pi*Fs),nip);
% for ii = 1:length(puls); DATA(1,1,ii) = h4c(1i*puls(ii))-D; end
% [Hloe,info] = mdspack.loewner(puls, DATA, struct('stable',true));
% Hloe.D      = D;
% mdspack.nyquist(Fc,puls1,'-')

% step 2: spectral zeros
tol_sz          = 1e-10;
[sz_R,sz_la]    = lf.spectral_zeros(Fcd);
sz_la_pos       = sz_la(real(sz_la)>tol_sz);
sz_R_pos        = sz_R(:,real(sz_la)>tol_sz);
h_fcd           = @(s) Fcd.c*((Fcd.e*s-Fcd.a)\Fcd.b) + Fcd.d;
%figure, plot(real(sz_la),imag(sz_la),'*')

% step 3: spectral zeros Loewner interpolation
clear la R W
for ii = 1:numel(sz_la_pos)
    la(ii)      = sz_la_pos(ii);
    R(:,ii)     = sz_R_pos(end-1+1:end,ii)/norm(sz_R_pos(end-1+1:end,ii));
    W(1,1,ii)   = h_fcd(sz_la_pos(ii));
end
%opt.D           = D;
opt.target      = [];
[~,info2]       = lf.loewner_tng(la,-conj(la),W,-conj(W),R,R',opt);
% Hloe     = info2.Hr;
% %Hloe        = c2d(info2.Hr,Ts,'tustin');
% h_loe       = @(z) Hloe.c*((Hloe.e*z-Hloe.a)\Hloe.b) + Hloe.d;
% 
% %tol_hermite = 1e-10;
% %[isstable(info2.Hr) isPassive(info2.Hr) norm(info2.LL-conj(info2.LL.'))<tol_hermite]
% % step 4: normalized passive
% [T,flag]    = chol(info2.LL); % R'*R = E (E > 0)
% if flag == 0
%     invT    = T\eye(length(T));
% else
%     warning('Cholesky decomposition issue')
%     invT    = eye(length(info2.Hr.A));
% end
% Hr          = ss(-invT'*info2.Hr.A*invT,-invT'*info2.Hr.B,info2.Hr.C*invT,info2.Hr.D);
% Hr        = c2d(Hr,1,'tustin');
% rr          = length(Hr.A);
% In          = eye(rr);
% hrp         = @(s) Hr.C*((s*In-Hr.A)\Hr.B)+Hr.D;
% 
% %[hloeph,info_loeph] = lf.passive2ph(Hr);


%%% Show results
figure, 
subplot(2,3,[1 2 4 5]), hold on, grid on
contour(X,Y,log10(abs(Zr3_loe)),lev_contour,'LineWidth',1,'DisplayName','$\mathbf{h}_3$')%,'ShowText','on')
plot(real(data.E),imag(data.E),'.','Color',[1 1 1]*.4,'MarkerSize',mw,'DisplayName',['(E) ' num2str(min(data.bnd),'%+2.0f')])
plot(real(data.F),imag(data.F),'k.','MarkerSize',mw,'DisplayName',['(F) ' num2str(max(data.bnd),'%+2.0f')])
plot(real(h4poles),imag(h4poles),'o','Color',col1,'MarkerFaceColor',col1,'DisplayName','(Z4) $p(\mathbf{h}_4)$')
plot(real(h4zeros),imag(h4zeros),'o','Color',col2,'DisplayName','(Z4) $z(\mathbf{h}_4)$')
plot(real(h3poles),imag(h3poles),'ro','MarkerFaceColor','r','DisplayName','(Z3) $p(\mathbf{h}_3)$')
plot(real(h3zeros),imag(h3zeros),'bo','MarkerFaceColor','b','DisplayName','(Z3) $z(\mathbf{h}_3)$')
plot(real(eigF),imag(eigF),'s','MarkerSize',mw,'DisplayName','IIR filter $p(\mathbf F)$')
plot(cos(2*pi*linspace(0,1,1e3)),sin(2*pi*linspace(0,1,1e3)),'k-','LineWidth',1,'DisplayName','$\mathcal C(0,1)$')
colormap winter
axis equal, set(gca,'Xlim',data.Xlim,'YLim',data.Ylim)
ylabel('Imag(.)'), xlabel('Real(.)')
title(['$r=' num2str(robj) '$, $\sigma_r=$ ' num2str(abs(hsig)) ' in ' num2str(timeLOE) 's' ])
legend('show','Location','eastoutside');
drawnow
%
w_val   = angle(pts);
w_chk   = linspace(-pi,pi,1e3);
z_chk   = exp(1i*w_chk);
for i = 1:numel(z_chk); G(i) = h4(z_chk(i)); end
for i = 1:numel(z_chk); Gloe(i) = hloeph(z_chk(i)); end
subplot(2,3,3), hold on, grid on
plot(w_val,20*log10(abs(val.')),'.','MarkerSize',20,'DisplayName','Data ($E \cup F$)'), hold on
plot(w_chk,20*log10(abs(G.'))  ,'-','LineWidth',2,'DisplayName','Filter (stable)')
plot(w_chk,20*log10(abs(Gloe.')),'--','LineWidth',2,'DisplayName','pH')
legend('show','Location','best')
ylabel('Gain [dB]'), xlabel('Frequency'), 
title(['Filter $\mathbf F$ gain $r=' num2str(robj) '$' ])
subplot(2,3,6), hold on, grid on
plot(w_val,angle(val),'.','MarkerSize',20,'DisplayName','Data ($E \cup F$)'), hold on
plot(w_chk,angle(G)  ,'-','LineWidth',2,'DisplayName','Filter (stable)')
%plot(w_chk,angle(Gloe)  ,'-','LineWidth',2,'DisplayName','pH')
title(['Filter $\mathbf F$ phase $r=' num2str(robj) '$' ])
ylabel('Phase [rad]'), xlabel('Frequency'), 
%
sgtitle([methodName ' filter design'],'FontSize',20)
drawnow
if SAVEIT, zol.figSavePDF('tex_pdf/figures/filter/iir_filter-pH'), end
%%
%%% Time-domain signal filtering
% >> generate signal and filter
Ns      = 2^14;
opt     = struct('f_band',[0 .9*Fs/2]);
[u,t]   = mdspack.signals('CHRP',Ns,Ts,opt);
y       = lsim(Fz,u);
% >> FFT
nt  = length(t);
f   = linspace(0,1,nt).'/Ts;
U   = fft(u)/nt;
Y   = fft(y)/nt;
% 
figure, 
subplot(211),hold on, axis tight
plot(t,u,'DisplayName','Original signal')
plot(t,y,'DisplayName','Filtered signal')
xlabel('Time [s]'),
%plot(t,u-y)
subplot(212),hold on, axis tight
plot(f,abs(U),'DisplayName','Original signal')
plot(f,abs(Y),'DisplayName','Filtered signal')
legend('show')
xlabel('Frequency [Hz]'), ylabel('Amplitude')
%plot(f,abs(U-Y))
sgtitle([methodName ' filter applied to chirp signal'],'FontSize',20)
if SAVEIT, zol.figSavePDF('tex_pdf/figures/filter/iir_filter-pH_td'), end

