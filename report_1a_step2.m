clearvars; close all; clc; format shorte %shortEng
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN: TO ADAPT BY USER
% Add Zolotarev Loewner package
addpath('/Users/charles/Documents/GIT/zolotarev')
% Add AAA package
addpath('/Users/charles/Documents/GIT/_others/chebfun')
% Directory name
fileDir = 'AlexReis_PCWIN64perso_MR2017a';
%%% END: TO ADAPT BY USER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Create folders
fileData    = ['tex_pdf' filesep 'data' filesep fileDir filesep];
fileFig     = [fileData 'figures' filesep];
fileTab     = [fileData 'tables' filesep];
allFiles    = zol.getFileData(fileData);
if ~exist(fileData); mkdir(fileData); end
if ~exist(fileFig);  mkdir(fileFig);  end
if ~exist(fileTab);  mkdir(fileTab);  end

%%% Style & slack variables
set(groot,'DefaultFigurePosition', [200 150 1000 600]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',6)
set(groot,'defaultaxesfontsize',20)
set(groot,'defaultAxesTickLabelInterpreter','latex');  
list_facto  = fieldnames(get(groot,'factory'));index_interpreter = find(contains(list_facto,'Interpreter'));for i = 1:length(index_interpreter); set(groot, strrep(list_facto{index_interpreter(i)},'factory','default'),'latex'); end
mw          = 20;
lw          = 4;
col         = parula(10);
col(1,:)    = [0 0 0];
col(2,:)    = [1 0 0];
mark        = {'s' 'o' '<' '>' '^' 'v'};
fileDirName = strrep(fileDir,'_',' ');

%%% Start loop
syms z
for i = 1:numel(allFiles)
    %%% Load file
    [fileData allFiles{i}]
    load([fileData allFiles{i}])
    N4 = []; D4 = []; ZE4 = []; PO4 = []; 
    %%% NUM/DEN ZER/POL EVAL
    close all
    for j = 1:numel(Z4)
        if isempty(NUM{j})
            [n4,d4] = zol.get_numden(Z4{j}); 
        else
            n4      = NUM{j};
            d4      = DEN{j};
        end
        N4  = [N4 [zeros(robj-numel(n4),1); n4]];
        D4  = [D4 [zeros(robj-numel(d4),1); d4]];
        ZE4 = [ZE4 [NaN*ones(robj-numel(ZER{j}),1); ZER{j}]];
        PO4 = [PO4 [NaN*ones(robj-numel(POL{j}),1); POL{j}]];
        for jj = 1:length(data.E)
            EVAL_E(jj,j) = Z4{j}(data.E(jj));
            EVAL_F(jj,j) = Z4{j}(data.F(jj));
        end
        %%%
        figure(1), 
        subplot(121), hold on, grid on
        plot(real(ZER{j}),imag(ZER{j}),mark{j},'Color',col(j,:),'MarkerSize',mw,'LineWidth',lw,'DisplayName',name{j})
        xlabel('Real'); ylabel('Imag.'); title('Zeros')
        legend('show','Location','best');
        subplot(122), hold on, grid on
        plot(real(POL{j}),imag(POL{j}),mark{j},'Color',col(j,:),'MarkerSize',mw,'LineWidth',lw,'DisplayName',name{j})
        xlabel('Real'); ylabel('Imag.'); title('Poles')
        sgtitle(['Approximation $r=' num2str(robj) '$'],'Fontsize',20)
    end
    %%% Figures
    figure(1)
    zol.figSavePDF([fileFig 'pz_r' num2str(robj)]), drawnow, pause(.5)
    %
    zer_max = 0; pol_max = 0;
    for k = 2:length(ZER)
        zer_max = max(zer_max,max(abs(real(ZER{k}(1:end-1)))));
        pol_max = max(pol_max,max(abs(real(POL{k}(1:end-1)))));
    end
    ZER_XLIM = [-1 1]*zer_max*1.5;%max(real(ZER{2}))*1e3;
    POL_XLIM = [-1 1]*pol_max*1.5;%max(real(ZER{2}))*1e3;
    figure(1)
    subplot(121), hold on,
    xlim(ZER_XLIM); ylim([-1 1]*max(imag(ZER{1}))*1.5); 
    subplot(122), hold on,
    xlim(POL_XLIM); ylim([-1 1]*max(imag(POL{1}))*1.5); 
    zol.figSavePDF([fileFig 'pz_r' num2str(robj) '_zoom']), drawnow, pause(.5)
    %
    figure, 
    for kk = 1:size(EVAL_E,2)
        subplot(121); hold on, axis equal, grid on, 
        plot(real(EVAL_E(:,kk)),imag(EVAL_E(:,kk)),'-','Color',col(kk,:))
        title('$E$')
        xlabel('Real'); ylabel('Imag.');
        subplot(122); hold on, axis equal, grid on, 
        plot(real(EVAL_F(:,kk)),imag(EVAL_F(:,kk)),'-','Color',col(kk,:))
        title('$F$')
        xlabel('Real'); ylabel('Imag.');
    end
    legend(name,'Location','South')
    sgtitle('Case 1a functions evaluation','Fontsize',20)
    zol.figSavePDF([fileFig 'eval_r' num2str(robj)]), drawnow, pause(.5)
    %%% LATEX
    % Num/Den
    Bnum        = (z.^(numel(n4)-1:-1:0)).';
    Bden        = (z.^(numel(d4)-1:-1:0)).';
    filename    = [fileTab '1a_nd_r' num2str(robj) '.tex'];
    str_coeff   = zol.latex_tab_nd(Bnum,Bden,N4,D4,name);
    str_coeff   = strrep(str_coeff,'\','\\');
    fileID      = fopen(filename, 'w');
    fprintf(fileID, '\\begin{table}[H] \\tiny $$ \n');
    fprintf(fileID, str_coeff);
    fprintf(fileID, [' \n $$\\normalsize \\caption{Case \\texttt{1a} (\\texttt{' fileDirName '}), $r=' num2str(robj) '$, Z4: numerator (first lines) and denominator (last lines) coefficients.} \\label{tab:sym-1a-r' num2str(robj)  '-' fileDir '} \\end{table}']);
    fclose(fileID);
    % Poles
    filename    = [fileTab '1a_pol_r' num2str(robj) '.tex'];
    fileID      = fopen(filename, 'w');
    str_pz      = zol.latex_tab_pz(PO4,name);
    fprintf(fileID, '\\begin{table}[H] \\tiny $$ \n');
    fprintf(fileID, str_pz);
    fprintf(fileID, ['\n $$\\normalsize \\caption{Case \\texttt{1a} (\\texttt{' fileDirName '}), $r=' num2str(robj) '$, Z4: poles.} \\label{tab:pol-1a-r' num2str(robj) '-' fileDir '} \\end{table}']);
    fclose(fileID);
    % Zeros
    filename    = [fileTab '1a_zer_r' num2str(robj) '.tex'];
    fileID      = fopen(filename, 'w');
    str_pz      = zol.latex_tab_pz(ZE4,name);
    fprintf(fileID, '\\begin{table}[H] \\tiny $$ \n');
    fprintf(fileID, str_pz);
    fprintf(fileID, ['\n $$\\normalsize \\caption{Case \\texttt{1a} (\\texttt{' fileDirName '}), $r=' num2str(robj) '$, Z4: zeros.} \\label{tab:zer-1a-r' num2str(robj)  '-' fileDir '} \\end{table}']);
    fclose(fileID);
end
license('inuse')
