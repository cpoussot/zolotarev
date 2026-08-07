clearvars; close all; clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN: TO ADAPT BY USER
% Add Zolotarev Loewner package
addpath('/Users/charles/Documents/GIT/zolotarev')
% Add AAA package
addpath('/Users/charles/Documents/GIT/_others/chebfun')
%%% Cases
spaceCAS    = {'1a' '1b' '1c' '1d' '1e' '1f' ...
               '2a' '2b' '2c' '2d' ...
               '3a' '3b' '3c' '3d' ...
               '7' 'spiral1' 'pm2'};
SHOW_ADVANC = false;
%%% END: TO ADAPT BY USER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Plot properties
set(groot,'DefaultFigurePosition', [200 150 1000 600]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',4)
set(groot,'defaultaxesfontsize',18)
set(groot,'defaultAxesTickLabelInterpreter','latex');  
list_facto  = fieldnames(get(groot,'factory'));index_interpreter = find(contains(list_facto,'Interpreter'));for i = 1:length(index_interpreter); set(groot, strrep(list_facto{index_interpreter(i)},'factory','default'),'latex'); end

%%% Loop over all cases
hsig_ = Inf*ones(numel(spaceCAS),100);
rsig_ = hsig_;
%figure
timeLoewner = 0;
timeAAA     = 0;
for j = 1:10%numel(spaceCAS)
    CAS = spaceCAS{j}
    %%% Define Zolotarev topology
    [pts,val,data]  = zol.example(CAS);
    %%% Estimate approximation bounds (thanks to Loewner)
    [la,mu,W,V]     = zol.example2data(pts,val,data);
    opt.target      = 1e-16;
    [h4,info]       = zol.loewner(la,mu,W,V,opt);
    rmax(j)         = info.r;
    %%% Now loop for all rational orders
    for i = 1:floor(rmax(j)*1.25)
        robj = i;
        %%% Loewner
        % (Z3-Z4)
        tic
        [la,mu,W,V]     = zol.example2data(pts,val,data);
        opt.target      = robj;
        [h4,info]       = zol.loewner(la,mu,W,V,opt);
        [h3,hp,hsig]    = zol.pb4_to_pb3(h4,pts,val);
        robj            = info.r;
        hsig_(j,i)      = hsig;
        timeLoewner     = timeLoewner + toc;

        %%% AAA 
        % (Z3-Z4)
        tic
        r4              = aaa(val,pts,"degree",robj,'sign',1,'damping',.95,'lawson',200);
        [r3,rp,rsig]    = zol.pb4_to_pb3(r4,pts,val);
        rsig_(j,i)      = rsig;
        timeAAA         = timeAAA + toc;
        % Show advancement
        if SHOW_ADVANC
            subplot(121), hold on, grid on, axis tight
            imagesc(real(log10(abs(hsig_.'))))
            c = colormap('parula'); c(end,:) = [1 1 1]; colormap(c)
            colorbar; c2 = clim;
            stairs((1:j)-.5,rmax,'r:','LineWidth',3)
            set(gca,'TickLabelInterpreter','latex','FontSize',16)
            title('\bf{LF}','Interpreter','latex','FontSize',18)
            xlabel('Case','Interpreter','latex','FontSize',18)
            ylabel('Degree $r$','Interpreter','latex','FontSize',18)
            xticks(1:numel(spaceCAS)), xticklabels(spaceCAS)
            %
            subplot(122), hold on, grid on, axis tight
            imagesc(real(log10(abs(rsig_.'))))
            colorbar; c1 = clim; c3 = [min([c1 c2]), max([c1 c2])]; clim(c3)
            set(gca,'TickLabelInterpreter','latex','FontSize',16)
            title('\bf{AAA}','Interpreter','latex','FontSize',18)
            xlabel('Case','Interpreter','latex','FontSize',18)
            ylabel('Degree $r$','Interpreter','latex','FontSize',18)
            xticks(1:numel(spaceCAS)), xticklabels(spaceCAS)
            %
            subplot(121), colorbar off
            sgtitle('$\log_{10}(\sigma_r)$ for different cases and orders','FontSize',24)
            %
            drawnow
        end    
    end
    [timeLoewner timeAAA]
end

license('inuse')
%%
%%% Plot
figure
c   = colormap('turbo'); c(end,:) = [1 1 1]; colormap(c)
subplot(121), hold on, grid on, axis tight
imagesc(real(log10(abs(hsig_.'))))
colorbar; c2 = clim;
stairs((1:j+1)-.5,[rmax rmax(end)],':','Color','r','LineWidth',3,'DisplayName','Loewner aut. order')
set(gca,'TickLabelInterpreter','latex','FontSize',16)
title(['\bf{LF (in $' num2str(timeLoewner/60) '$ min)}'],'Interpreter','latex','FontSize',18)
xlabel('\bf{Case}','Interpreter','latex','FontSize',18)
ylabel('\bf{Degree $r$}','Interpreter','latex','FontSize',18)
xticks(1:numel(spaceCAS)), xticklabels(spaceCAS)
xlim([.5 j+.5])
legend('show')
%
subplot(122), hold on, grid on, axis tight
imagesc(real(log10(abs(rsig_.'))))
colorbar; c1 = clim; c3 = [min([c1 c2]), max([c1 c2])]; clim(c3)
stairs((1:j+1)-.5,[rmax rmax(end)],':','Color','r','LineWidth',3,'DisplayName','Loewner aut. order')
set(gca,'TickLabelInterpreter','latex','FontSize',16)
title(['\bf{AAA (in $' num2str(timeAAA/60) '$ min)}'],'Interpreter','latex','FontSize',18)
xlabel('\bf{Case}','Interpreter','latex','FontSize',18)
ylabel('\bf{Degree $r$}','Interpreter','latex','FontSize',18)
xticks(1:numel(spaceCAS)), xticklabels(spaceCAS)
xlim([.5 j+.5])
legend('show')
%
subplot(121), colorbar off
sgtitle('\bf{$\log_{10}(\sigma_r)$ for different cases and orders}','FontSize',24)
drawnow
zol.figSavePDF('tex_pdf/figures/all/case_oder')