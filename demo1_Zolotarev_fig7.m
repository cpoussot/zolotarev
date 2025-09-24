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
mw = 15;
warning('off')
spaceCAS    = {'1a' '1b' '1c' '1d' '1e' '1f' ...
               '2a' '2b' '2c' '2d' ...
               '3a' '3b' '3c' '3d' ...
               '7' 'spiral1' 'pm' 'pm2'};
lgn         = {'Loewner','AAA','AAA  \texttt{"sign",1}','AAA  \texttt{"sign",1,"damping",.95,"lawson",200}'};
%lgn{length(lgn)+1} = ['(E) ' num2str(-1,'%+2.0f')];
%lgn{length(lgn)+1} = ['(F) ' num2str(+1,'%+2.0f')];
for j = 1%:numel(spaceCAS)%l-1:numel(spaceCAS)
    close all
    clear hsig_ rsig1_ rsig2_ rsig3_ 
    clear timeLOE_ timeAAA1_ timeAAA2_ timeAAA3_ 
    CAS = spaceCAS{j}
    %%% Define Zolotarev topology
    [pts,val,data]  = zol.example(CAS);
    %%% Estimate bound
    [la,mu,W,V]     = zol.example2data(pts,val,data);
    opt.target      = 1e-16;
    [h4,info]       = zol.loewner(la,mu,W,V,opt);
    %
    figure, 
    for i = 1:info.r
        robj = i;
        %%% Loewner
        % (Z3-Z4)
        tic
        [la,mu,W,V]     = zol.example2data(pts,val,data);
        opt.target      = robj;
        [h4,info]       = zol.loewner(la,mu,W,V,opt);
        [h3,hp,hsig]    = zol.pb4_to_pb3(h4,pts,val);
        timeLOE         = toc;
        robj            = info.r;
    
        %%% AAA
        % (Z3-Z4)
        tic
        [r4,rpoles,~,rzeros,zj,fj,wj] = aaa(val,pts,"degree",robj);
        [r3,rp,rsig1]   = zol.pb4_to_pb3(r4,pts,val);
        timeAAA1        = toc;
        % (Z3-Z4)
        tic
        [r4,rpoles,~,rzeros,zj,fj,wj] = aaa(val,pts,"degree",robj,'sign',1);
        [r3,rp,rsig2]   = zol.pb4_to_pb3(r4,pts,val);
        timeAAA2        = toc;
        % (Z3-Z4)
        tic
        [r4,rpoles,~,rzeros,zj,fj,wj] = aaa(val,pts,"degree",robj,'sign',1,'damping',.95,'lawson',200);
        [r3,rp,rsig3]   = zol.pb4_to_pb3(r4,pts,val);
        timeAAA3        = toc;
    
        timeLOE_(i)     = timeLOE;
        timeAAA1_(i)    = timeAAA1;
        timeAAA2_(i)    = timeAAA2;
        timeAAA3_(i)    = timeAAA3;
        hsig_(i)        = hsig;
        rsig1_(i)       = rsig1;
        rsig2_(i)       = rsig2;
        rsig3_(i)       = rsig3;
        %
        clf
        subplot(2,3,[1 2]), hold on, grid on, axis tight
        plot(1:i,abs(hsig_),'-','LineWidth',3,'DisplayName',lgn{1})
        plot(1:i,abs(rsig1_),'--','LineWidth',3,'DisplayName',lgn{2})
        plot(1:i,abs(rsig2_),'-.','LineWidth',3,'DisplayName',lgn{3})
        plot(1:i,abs(rsig3_),':','LineWidth',3,'DisplayName',lgn{4})
        set(gca,'YScale','log')
        %legend(lgn,'Interpreter','latex','Location','SouthWest','FontSize',16)
        set(gca,'TickLabelInterpreter','latex','FontSize',16)
        title('\bf{Methods computed Zolotarev ratio $\sigma_r$}','Interpreter','latex','FontSize',18)
        ylabel('$\sigma_r$','Interpreter','latex','FontSize',18)
        xlabel('Degree $r$','Interpreter','latex','FontSize',18)
        subplot(2,3,[4 5]), hold on, grid on, axis tight
        plot(1:i,timeLOE_,'-','LineWidth',3)
        plot(1:i,timeAAA1_,'--','LineWidth',3)
        plot(1:i,timeAAA2_,'-.','LineWidth',3)
        plot(1:i,timeAAA3_,':','LineWidth',3)
        set(gca,'YScale','log')
        set(gca,'TickLabelInterpreter','latex','FontSize',16)
        ylabel('Time [s]','Interpreter','latex','FontSize',18)
        xlabel('Degree $r$','Interpreter','latex','FontSize',18)
        title('\bf{Methods computational time}','Interpreter','latex','FontSize',20)
        subplot(2,3,3), hold on, grid on, axis equal
        plot(real(data.E),imag(data.E),'.','Color',[1 1 1]*.4,'MarkerSize',mw,'DisplayName',['(E) ' num2str(min(data.bnd),'%+2.0f')])
        plot(real(data.F),imag(data.F),'k.','MarkerSize',mw,'DisplayName',['(F) ' num2str(max(data.bnd),'%+2.0f')])
        set(gca,'Xlim',data.Xlim,'YLim',data.Ylim)
        set(gca,'TickLabelInterpreter','latex','FontSize',16)
        ylabel('Imag(.)','Interpreter','latex','FontSize',16)
        xlabel('Real(.)','Interpreter','latex','FontSize',16)
        title('\bf{Geometry}','Interpreter','latex','FontSize',18)
        subplot(236), hold on
        plot(0,0,'-')
        plot(0,0,'--')
        plot(0,0,'-.')
        plot(0,0,':')
        %plot(0,0,'.','Color',[1 1 1]*.4,'MarkerSize',mw)
        %plot(0,0,'k.','MarkerSize',mw)
        h4 = gca;
        h4.XAxis.Visible = 'off';
        h4.YAxis.Visible = 'off';
        %set(gca,'xtick',[],'ytick',[])
        legend(lgn,'Interpreter','latex','Location','South','FontSize',12)
        % add a bit space to the figure
        drawnow
    end
    % fig = gcf;
    % fig.Position(3) = fig.Position(3) + 250;
    % Lgnd = legend('show');
    % Lgnd.Position(1) = 0.015;
    % Lgnd.Position(2) = 0.4;
    % drawnow
    lf.figSavePDF(['figures/time_accuracy/cas_' CAS '_zol_time'],.6)
end

license('inuse')
