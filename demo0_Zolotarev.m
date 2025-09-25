clearvars; close all; clc
set(groot,'DefaultFigurePosition', [200 150 1000 600]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',4)
set(groot,'defaultaxesfontsize',18)
set(groot,'defaultAxesTickLabelInterpreter','latex');  
list_factory = fieldnames(get(groot,'factory'));index_interpreter = find(contains(list_factory,'Interpreter'));for i = 1:length(index_interpreter); set(groot, strrep(list_factory{index_interpreter(i)},'factory','default'),'latex'); end
%%% AAA package
addpath('/Users/charles/Documents/GIT/chebfun')
%%% Chose case, order
CAS     = '1a' % /!\ '1a' and '1b' use "Symbolic Toolbox" 
robj0   = 1e-14;
mw      = 15; % marker width
kk      = 1; % for plot
%%% Chose options for chefun AAA
AAAparam    = {""; ... 
               ",'sign',1,'damping',.95"; ...
               ",'sign',1,'damping',.95,'lawson',200"};
K           = ceil((1+numel(AAAparam))/2);
lev_contour = -30:1:0;

%%% Define Zolotarev topology
[pts,val,data]  = zol.example(CAS);
xx              = linspace(data.Xlim(1),data.Xlim(2),101);
yy              = linspace(data.Ylim(1),data.Ylim(2),103);
[X,Y]           = meshgrid(xx,yy);

%%% Loewner approximation
% >> (Z4) rational approximation
tic
methodName      = '\textbf{LF}';
[la,mu,W,V]     = zol.example2data(pts,val,data);
opt             = [];
opt.target      = robj0;
opt.D           = 0;
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
Zr3_loe         = zeros(numel(yy),numel(xx));
for i = 1:numel(xx)
    for j = 1:numel(yy)
        Zr3_loe(j,i) = h3(xx(i)+1i*yy(j));
    end
end 
%
figure(1)
subplot(K,2,kk), hold on, grid on
contour(X,Y,log10(abs(Zr3_loe)),lev_contour,'LineWidth',1,'DisplayName','$\mathbf{h}_3$')%,'ShowText','on')
plot(real(data.E),imag(data.E),'.','Color',[1 1 1]*.4,'MarkerSize',mw,'DisplayName',['(E) ' num2str(min(data.bnd),'%+2.0f')])
plot(real(data.F),imag(data.F),'k.','MarkerSize',mw,'DisplayName',['(F) ' num2str(max(data.bnd),'%+2.0f')])
plot(real(h4poles),imag(h4poles),'ko','MarkerFaceColor','k','DisplayName','(Z4) $p(\mathbf{h}_4)$')
plot(real(h4zeros),imag(h4zeros),'ko','DisplayName','(Z4) $z(\mathbf{h}_4)$')
plot(real(h3poles),imag(h3poles),'ro','MarkerFaceColor','r','DisplayName','(Z3) $p(\mathbf{h}_3)$')
plot(real(h3zeros),imag(h3zeros),'bo','MarkerFaceColor','b','DisplayName','(Z3) $z(\mathbf{h}_3)$')
colormap winter
axis equal, set(gca,'Xlim',data.Xlim,'YLim',data.Ylim)
ylabel('Imag(.)'), xlabel('Real(.)')
title({methodName; ['$r=$' num2str(robj) ', $\sigma_r=$ ' num2str(abs(hsig)) ' in ' num2str(timeLOE) 's' ]})
Lgnd = legend('show');
Lgnd.Position(1) = 0.015;
Lgnd.Position(2) = 0.4;
drawnow

%%% AAA approximation for different options
for ii = 1:numel(AAAparam)
    kk              = kk + 1;
    AAAparam_i      = AAAparam{ii};
    % >> (Z4) rational approximation + poles and zeros
    tic
    eval(['[r4,r4poles,~,r4zeros,zj,fj,wj] = aaa(val,pts,"degree",robj' AAAparam_i{1} ');'])
    timeAAA         = toc;
    % >> (Z3) rational approximation, from Z4->Z3
    [r3,rp,rsig]    = zol.pb4_to_pb3(r4,pts,val);
    % >> (Z3) poles and zeros
    [~,~,r3poles]   = prz(zj,fj-rp,wj);
    [~,~,r3zeros]   = prz(zj,fj+rp,wj);
    % >> (Z3) evaluate
    Zr3_aaa         = zeros(numel(yy),numel(xx));
    for i = 1:numel(xx)
        for j = 1:numel(yy)
            Zr3_aaa(j,i) = r3(xx(i)+1i*yy(j));
        end
    end
    %
    figure(1), 
    subplot(K,2,kk), hold on, grid on
    contour(X,Y,log10(abs(Zr3_aaa)),lev_contour,'LineWidth',1)%,'ShowText','on')
    plot(real(data.E),imag(data.E),'.','Color',[1 1 1]*.4,'MarkerSize',mw)
    plot(real(data.F),imag(data.F),'k.','MarkerSize',mw)
    plot(real(r4poles),imag(r4poles),'ko','MarkerFaceColor','k')
    plot(real(r4zeros),imag(r4zeros),'ko')
    plot(real(r3poles),imag(r3poles),'ro','MarkerFaceColor','r')
    plot(real(r3zeros),imag(r3zeros),'bo','MarkerFaceColor','b')
    colormap winter
    axis equal, set(gca,'Xlim',data.Xlim,'YLim',data.Ylim)
    ylabel('Imag(.)'), xlabel('Real(.)')
    title({['\textbf{AAA} \texttt{opt=' AAAparam_i{1}(2:end) '}']; ...
           ['$r=$' num2str(robj) '(' num2str(length(h3poles)) '), $\sigma_r=$ ' num2str(abs(rsig)) ' in ' num2str(timeAAA) 's' ]})
    drawnow
end
license('inuse')
