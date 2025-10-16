clearvars; close all; clc
set(groot,'DefaultFigurePosition', [200 150 1000 600]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',4)
set(groot,'defaultaxesfontsize',18)
set(groot,'defaultAxesTickLabelInterpreter','latex');  
list_factory = fieldnames(get(groot,'factory'));index_interpreter = find(contains(list_factory,'Interpreter'));for i = 1:length(index_interpreter); set(groot, strrep(list_factory{index_interpreter(i)},'factory','default'),'latex'); end
%%% Zolotarev Loewner package
addpath('/Users/charles/Documents/GIT/zolotarev')
%%% AAA package
addpath('/Users/charles/Documents/GIT/chebfun')
%
CAS         = '1b'
robj        = 6;
Nplot       = 1;
interlace   = true;
AAAparam    = {""; ... 
               ",'sign',1,'damping',.95,'lawson',200"};
%%% Define Zolotarev topology
[pts,val,data]  = zol.example(CAS);
xx              = linspace(data.Xlim(1),data.Xlim(2),101);
yy              = linspace(data.Ylim(1),data.Ylim(2),103);
[X,Y]           = meshgrid(xx,yy);

%%% Optimal
syms z
n4_opt              = data.z4{robj}{1};
d4_opt              = data.z4{robj}{2};
h4_opt              = matlabFunction(poly2sym(n4_opt,z)/poly2sym(d4_opt,z));
[h3_opt,~,sig_opt]  = zol.pb4_to_pb3(h4_opt,pts,val);
deg_num             = length(n4_opt)-1;
deg_den             = length(d4_opt)-1;

%%% Loewner
[la,mu,W,V]         = zol.example2data(pts,val,data);
opt.target          = robj;
%opt.D               = 1;
[h4_loe,info]       = zol.loewner(la,mu,W,V,opt);
[n4_loe,d4_loe]     = tfdata(dss(info.Ar,info.Br,info.Cr,info.Dr,info.Er)); n4_loe = n4_loe{1}.'; d4_loe = d4_loe{1}.';
n4_loe              = [zeros(deg_num-length(n4_loe)+1,1); n4_loe];
d4_loe              = [zeros(deg_den-length(d4_loe)+1,1); d4_loe];
[h3_loe,~,sig_loe]  = zol.pb4_to_pb3(h4_loe,pts,val);
% Coefficients n/d
fprintf('*** Loewner Framework (LF) *** \n')
fprintf('Z4: numerator, optimal (left) LF (right)\n')
vpa([n4_opt n4_loe],3)
fprintf('Z4: denominator, optimal (left) LF (right)\n')
vpa([d4_opt d4_loe],3)
fprintf('************************ \n')
%
SIG = sig_loe;
NUM = n4_loe;
DEN = d4_loe;

%%% AAA 
ii = 1;
AAAparam_i  = AAAparam{ii};
eval(['[h4_aaa,h4_aaa_pol,h4_aaa_res,h4_aaa_zer,zj,fj,wj] = aaa(val,pts,"degree",robj' AAAparam_i{1} ');'])
[n4_aaa,d4_aaa,M]   = zol.get_numden(h4_aaa);
n4_aaa              = [zeros(deg_num-length(n4_aaa)+1,1); n4_aaa];
d4_aaa              = [zeros(deg_den-length(d4_aaa)+1,1); d4_aaa];
[h3_aaa,~,sig_aaa]  = zol.pb4_to_pb3(h4_aaa,pts,val);
% Coefficients n/d
fprintf('*** AAA %s ***\n',AAAparam{ii})
fprintf('Z4: numerator, optimal (left) AAA (right)\n')
vpa([n4_opt n4_aaa],3)
fprintf('Z4: denominator, optimal (left) AAA (right)\n')
vpa([d4_opt d4_aaa],3)
fprintf('************************ \n')

SIG = [SIG sig_aaa];
NUM = [NUM n4_aaa];
DEN = [DEN d4_aaa];

%%% AAA 
ii = 2;
AAAparam_i  = AAAparam{ii};
eval(['[h4_aaa,h4_aaa_pol,h4_aaa_res,h4_aaa_zer,zj,fj,wj] = aaa(val,pts,"degree",robj' AAAparam_i{1} ');'])
[n4_aaa,d4_aaa,M]   = zol.get_numden(h4_aaa);
n4_aaa              = [zeros(deg_num-length(n4_aaa)+1,1); n4_aaa];
d4_aaa              = [zeros(deg_den-length(d4_aaa)+1,1); d4_aaa];
[h3_aaa,~,sig_aaa]  = zol.pb4_to_pb3(h4_aaa,pts,val);
% Coefficients n/d
fprintf('*** AAA %s ***\n',AAAparam{ii})
fprintf('Z4: numerator, optimal (left) AAA (right)\n')
vpa([n4_opt n4_aaa],3)
fprintf('Z4: denominator, optimal (left) AAA (right)\n')
vpa([d4_opt d4_aaa],3)
fprintf('************************ \n')

SIG = [SIG sig_aaa];
NUM = [NUM n4_aaa];
DEN = [DEN d4_aaa];



%%% Plot
%
Bnum    = (z.^(deg_num:-1:0)).';
Bden    = (z.^(deg_den:-1:0)).';
%
g4_opt  = matlabFunction(sum(n4_opt.*Bnum)/sum(d4_opt.*Bden));
g4_loe  = matlabFunction(sum(NUM(:,1).*Bnum)/sum(DEN(:,1).*Bden));
g4_aaa  = matlabFunction(sum(NUM(:,2).*Bnum)/sum(DEN(:,2).*Bden));
g4_aaaL = matlabFunction(sum(NUM(:,3).*Bnum)/sum(DEN(:,3).*Bden));

switch CAS
    case '1a'
        xx = data.F;
    case '1b'
        xx = sort([data.E data.F]);
end

for i = 1:length(xx)
    g4_opt_(i)  = g4_opt(xx(i));
    g4_loe_(i)  = g4_loe(xx(i));
    g4_aaa_(i)  = g4_aaa(xx(i));
    g4_aaaL_(i) = g4_aaaL(xx(i));
end

col         = parula(5);
col(2,:)    = [1, 0, 0];

switch CAS
    case '1a'
        figure, hold on, axis equal, grid on
        plot(real(g4_opt_),imag(g4_opt_),'-','Color',col(1,:))
        plot(real(g4_loe_),imag(g4_loe_),'-','Color',col(2,:))
        plot(real(g4_aaa_),imag(g4_aaa_),'--s','Color',col(3,:))
        plot(real(g4_aaaL_),imag(g4_aaaL_),'-.','Color',col(4,:))
        xlabel('Real'); ylabel('Imag.');
        title(['Case ' CAS ': functions evaluation'])
        legend({'Optimal','LF',strcat('AAA ',AAAparam{1}),strcat('AAA ',AAAparam{2})},'location','best')
        %legend({'Optimal','LF','AAA','AAA-Lawson'},'location','best')
        xlim([min(real(g4_loe_)) max(real(g4_loe_))]); 
        ylim([min(imag(g4_loe_)) max(imag(g4_loe_))]); 
    case '1b'
        figure
        subplot(211), hold on, axis tight, grid on
        plot(xx,g4_opt_,'-','Color',col(1,:))
        plot(xx,g4_loe_,'-','Color',col(2,:))
        plot(xx,g4_aaa_,'--s','Color',col(3,:))
        plot(xx,g4_aaaL_,'-.','Color',col(4,:))
        xlabel('$x$'); xlabel('$\mathbf{f}(x)$');
        title(['Case ' CAS ': functions evaluation'])
        legend({'Optimal','LF',strcat('AAA ',AAAparam{1}),strcat('AAA ',AAAparam{2})},'location','best')
        %legend({'Optimal','LF','AAA','AAA-Lawson'},'location','best')
        subplot(212), hold on, axis tight, grid on
        plot(xx,abs(g4_loe_-g4_opt_),'-','Color',col(2,:))
        plot(xx,abs(g4_aaa_-g4_opt_),'-s','Color',col(3,:))
        plot(xx,abs(g4_aaaL_-g4_opt_),'-.','Color',col(4,:))
        set(gca,'YScale','log')
        xlabel('$x$'); ylabel('$|\mathbf{f}(x)-\mathbf{f}^\star(x)|$');
        title(['Case ' CAS ': functions error evaluation'])
end
