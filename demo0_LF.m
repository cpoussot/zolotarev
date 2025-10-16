clearvars; close all; clc
set(groot,'DefaultFigurePosition', [200 150 1000 600]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',4)
set(groot,'defaultaxesfontsize',18)
set(groot,'defaultAxesTickLabelInterpreter','latex');  
list_factory = fieldnames(get(groot,'factory'));index_interpreter = find(contains(list_factory,'Interpreter'));for i = 1:length(index_interpreter); set(groot, strrep(list_factory{index_interpreter(i)},'factory','default'),'latex'); end
%%% Chose case, order
CAS         = '7';    % /!\ '1a' and '1b' use "Symbolic Toolbox" if available
robj0       = 1e-14;   % objective order (either integer > 1 or sigma threshold)

%%% Plot properties
mw          = 15;         % marker width
kk          = 1;          % for plot
lev_contour = -30:1:0;    % contour plot
col         = parula(10); % 
col1        = col(5,:);   % Z4 poles
col2        = col(9,:);   % Z4 zeros

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
figure, hold on
contour(X,Y,log10(abs(Zr3_loe)),lev_contour,'LineWidth',1,'DisplayName','$\mathbf{h}_3$')%,'ShowText','on')
plot(real(data.E),imag(data.E),'.','Color',[1 1 1]*.4,'MarkerSize',mw,'DisplayName',['(E) ' num2str(min(data.bnd),'%+2.0f')])
plot(real(data.F),imag(data.F),'k.','MarkerSize',mw,'DisplayName',['(F) ' num2str(max(data.bnd),'%+2.0f')])
plot(real(h4poles),imag(h4poles),'o','Color',col1,'MarkerFaceColor',col1,'DisplayName','(Z4) $p(\mathbf{h}_4)$')
plot(real(h4zeros),imag(h4zeros),'o','Color',col2,'DisplayName','(Z4) $z(\mathbf{h}_4)$')
plot(real(h3poles),imag(h3poles),'ro','MarkerFaceColor','r','DisplayName','(Z3) $p(\mathbf{h}_3)$')
plot(real(h3zeros),imag(h3zeros),'bo','MarkerFaceColor','b','DisplayName','(Z3) $z(\mathbf{h}_3)$')
colormap winter
axis equal, set(gca,'Xlim',data.Xlim,'YLim',data.Ylim)
ylabel('Imag(.)'), xlabel('Real(.)')
title({methodName; ['$r=' num2str(robj) '$, $\sigma_r=$ ' num2str(abs(hsig)) ' in ' num2str(timeLOE) 's' ]})
Lgnd = legend('show');
drawnow

license('inuse')
