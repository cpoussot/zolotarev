clearvars; close all; clc
set(groot,'DefaultFigurePosition', [200 150 1000 600]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',4)
set(groot,'defaultaxesfontsize',18)
set(groot,'defaultAxesTickLabelInterpreter','latex');  
list_factory = fieldnames(get(groot,'factory'));index_interpreter = find(contains(list_factory,'Interpreter'));for i = 1:length(index_interpreter); set(groot, strrep(list_factory{index_interpreter(i)},'factory','default'),'latex'); end
% %%% Zolotorev Loewner package
% addpath('/Users/charles/Documents/GIT/zolotorev')
%%% AAA package
addpath('/Users/charles/Documents/GIT/chebfun')
%
SAVEIT  = false; 
CAS     = '1d'
robj0   = 1e-14;
Nplot   = 1;
mw      = 15;
%
% CAS = '1a'; robj0=20;    SAVEIT = true;
% CAS = '1a'; robj0=1e-14; SAVEIT = true;
% CAS = '1a'; robj0=1e-15; SAVEIT = true;
% CAS = 'spiral1'; robj0=1e-15; SAVEIT = true;
% CAS = 'pm'; robj0=1e-16; SAVEIT = true;
% CAS = 'pm2'; robj0=1e-16; SAVEIT = true;
%
% AAAparam    = {""; ...
%                ",'sign',1,'damping',.95"; ...
%                ",'sign',1,'damping',.95,'lawson',10"; ...
%                ",'sign',1,'damping',.95,'lawson',200"; ...
%                ",'sign',1,'damping',.95,'lawson',400"};
AAAparam    = {""; ... 
               ",'sign',1,'damping',.95"; ...
               ",'sign',1,'damping',.95,'lawson',200"};
% AAAparam    = {"",",'sign',1,'damping',.95,'lawson',200"};
%K           = 1+ceil(numel(AAAparam)/2);
K           = ceil(numel(AAAparam)/2);
lev         = 30;
lev_contour = -lev:1:0;
kk          = 0;

%%% Define Zolotarev topology
[pts,val,data]  = zol.example(CAS);
xx              = linspace(data.Xlim(1)*Nplot,data.Xlim(2)*Nplot,101);
yy              = linspace(data.Ylim(1)*Nplot,data.Ylim(2)*Nplot,103);
[X,Y]           = meshgrid(xx,yy);

%%% Loewner
% >> (Z4) \hat h, computation + poles and zeros
for ii = 1%:2
    kk  = kk +1
    tic
    switch ii 
        case 1
            methodName  = '\textbf{LF}';
            [la,mu,W,V] = zol.example2data(pts,val,data);
            opt         = [];
            opt.target  = robj0;
            opt.D       = 0;
            [h4,info]   = zol.loewner(la,mu,W,V,opt);
            %[h4,info]   = zol.loewner_zol(val,pts,robj0,n1,n2,d);
        case 2
            methodName  = '\textbf{LF (with deriv.)}';
            %[la,mu,W,V] = zol.example2data(pts,val,data);
            % shift
            % tol = 1e-9;
            % W_shift = repmat([tol tol -tol -tol],1,length(la)/4);
            % V_shift = repmat([tol;tol;-tol;-tol],length(mu)/4,1);
            % W = W+W_shift;
            % V = V+V_shift;
            [h4,info]   = zol.loewner_deriv(val,pts,robj0,0,1e-19);
    end
    %
    test1 = (info.LL*info.LA - info.MU*info.LL);
    test2 = -(info.V*info.R - info.L*info.W);
    norm(test1-test2)
    % [norm(test1) norm(test2)]
    h4poles         = eig(info.Ar,info.Er);
    h4zeros         = eig([info.Ar info.Br;info.Cr 0],blkdiag(info.Er,0));
    timeLOE         = toc;
    robj            = info.r;
    % figure(ii), hold on
    % plot(la,'x')
    % plot(mu,'o')
    % drawnow
    % subplot(121), imagesc(log10(abs(info.LL))), colorbar
    % subplot(122), imagesc(log10(abs(info.SS))), colorbar
    
    % >> (Z3) h*, from Z4->Z3 computation + poles and zeros
    % h4 = @(s) sqrt(hsig)*(hp+info.C*((s*info.E-info.A)\info.B))/...
    %                      (hp-info.C*((s*info.E-info.A)\info.B));
    warning('off')
    [h3,hp,hsig]    = zol.pb4_to_pb3(h4,pts,val);
    h3poles         = eig([info.Ar info.Br;-info.Cr (hp)],blkdiag(info.Er,0));
    h3zeros         = eig([info.Ar info.Br; info.Cr (hp)],blkdiag(info.Er,0));
    % Symbolic to check coefficients
    %[h3num,h3den]   = zol.get_numden(h3,true);
    
    Zr3_loe = zeros(numel(yy),numel(xx));
    for i = 1:numel(xx)
        for j = 1:numel(yy)
            Zr3_loe(j,i) = h3(xx(i)+1i*yy(j));
        end
    end 
    
    figure(100)
    %subplot(K,2,1,'Color',[229 255 204]/255), hold on, grid on
    subplot(K,2,kk), hold on, grid on
    %surfc(X,Y,log10(abs(Zr3_loe)),'DisplayName','$\mathbf{h}_3$')%,'ShowText','on')
    contour(X,Y,log10(abs(Zr3_loe)),lev_contour,'LineWidth',1,'DisplayName','$\mathbf{h}_3$')%,'ShowText','on')
    plot(real(data.E),imag(data.E),'.','Color',[1 1 1]*.4,'MarkerSize',mw,'DisplayName',['(E) ' num2str(min(data.bnd),'%+2.0f')])
    plot(real(data.F),imag(data.F),'k.','MarkerSize',mw,'DisplayName',['(F) ' num2str(max(data.bnd),'%+2.0f')])
    plot(real(h4poles),imag(h4poles),'ko','MarkerFaceColor','k','DisplayName','(Z4) $p(\mathbf{h}_4)$')
    plot(real(h4zeros),imag(h4zeros),'ko','DisplayName','(Z4) $z(\mathbf{h}_4)$')
    plot(real(h3poles),imag(h3poles),'ro','MarkerFaceColor','r','DisplayName','(Z3) $p(\mathbf{h}_3)$')
    plot(real(h3zeros),imag(h3zeros),'bo','MarkerFaceColor','b','DisplayName','(Z3) $z(\mathbf{h}_3)$')
    %plot(real(info.lar),imag(info.lar),'s')
    %plot(real(info.mur),imag(info.mur),'x')
    colormap winter
    % legend({'$\mathbf{h}_3$', ... 
    %         ['(E) ' num2str(min(data.bnd),'%+2.0f')], ...
    %         ['(F) ' num2str(max(data.bnd),'%+2.0f')], ...
    %         '(Z4) $p(\mathbf{h}_4)$', ...
    %         '(Z4) $z(\mathbf{h}_4)$', ...
    %         '(Z3) $p(\mathbf{h}_3)$', ...
    %         '(Z3) $z(\mathbf{h}_3)$'}, ... 
    %         'Location','eastoutside')
    %         %'IP (right)','IP (left)'}
    axis equal, set(gca,'Xlim',data.Xlim*Nplot,'YLim',data.Ylim*Nplot)
    ylabel('Imag(.)'), xlabel('Real(.)')
    title({methodName; ['$r=$' num2str(robj) ', $\sigma_r=$ ' num2str(abs(hsig)) ' in ' num2str(timeLOE) 's' ]})
    % add a bit space to the figure
    %fig = gcf;
    %fig.Position(3) = fig.Position(3) + 250;
    Lgnd = legend('show');
    Lgnd.Position(1) = 0.015;
    Lgnd.Position(2) = 0.4;
    drawnow
end
%%
%%% AAA 
%kk = 1;
for ii = 1:numel(AAAparam)
    % >> (Z4) \hat r computation + poles and zeros
    kk          = kk + 1;
    AAAparam_i  = AAAparam{ii};
    tic
    eval(['[r4,r4poles,~,r4zeros,zj,fj,wj] = aaa(val,pts,"degree",robj' AAAparam_i{1} ');'])
    timeAAA     = toc;
    % >> (Z3) from Z4->Z3 + poles and zeros
    [r3,rp,rsig]    = zol.pb4_to_pb3(r4,pts,val);
    [~,~,r3poles]   = prz(zj,fj-rp,wj);
    [~,~,r3zeros]   = prz(zj,fj+rp,wj);
    % % Symbolic to check coefficients
    % [r3num,r3den]   = zol.get_numden(r3,true);
    %
    Zr3_aaa = zeros(numel(yy),numel(xx));
    for i = 1:numel(xx)
        for j = 1:numel(yy)
            Zr3_aaa(j,i) = r3(xx(i)+1i*yy(j));
        end
    end
    %
    figure(100), 
    subplot(K,2,kk), hold on, grid on
    contour(X,Y,log10(abs(Zr3_aaa)),lev_contour,'LineWidth',1)%,'ShowText','on')
    plot(real(data.E),imag(data.E),'.','Color',[1 1 1]*.4,'MarkerSize',mw)
    plot(real(data.F),imag(data.F),'k.','MarkerSize',mw)
    plot(real(r4poles),imag(r4poles),'ko','MarkerFaceColor','k')
    plot(real(r4zeros),imag(r4zeros),'ko')
    plot(real(r3poles),imag(r3poles),'ro','MarkerFaceColor','r')
    plot(real(r3zeros),imag(r3zeros),'bo','MarkerFaceColor','b')
    colormap winter
    % legend({'$\mathbf{r}_3$', ... 
    %         ['(E) ' num2str(min(data.bnd),'%+2.0f')], ...
    %         ['(F) ' num2str(max(data.bnd),'%+2.0f')], ...
    %         '(Z4) $p(\mathbf{r}_4)$', ...
    %         '(Z4) $z(\mathbf{r}_4)$', ...
    %         '(Z3) $p(\mathbf{r}_3)$', ...
    %         '(Z3) $z(\mathbf{r}_3)$'},'Location','eastoutside')
    axis equal, set(gca,'Xlim',data.Xlim*Nplot,'YLim',data.Ylim*Nplot)
    ylabel('Imag(.)'), xlabel('Real(.)')
    title({['\textbf{AAA} \texttt{opt=' AAAparam_i{1}(2:end) '}']; ...
           ['$r=$' num2str(robj) '(' num2str(length(h3poles)) '), $\sigma_r=$ ' num2str(abs(rsig)) ' in ' num2str(timeAAA) 's' ]})
    drawnow
end
if SAVEIT
    lf.figSavePDF(['figures/intro/cas_' CAS '_r' num2str(robj0) ],.6)
end
license('inuse')
