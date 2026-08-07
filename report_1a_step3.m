clearvars; close all; clc; format shorte
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN: TO ADAPT BY USER
% Add Zolotarev Loewner package
addpath('/Users/charles/Documents/GIT/zolotarev')
% Add AAA package
addpath('/Users/charles/Documents/GIT/_others/chebfun')
% Directory name list
fileDir = {'cpv_MACA64_MR2023b' ...
           'cpv_MACA64_MR2022b' ... 
           'cpv_PCWIN64_MR2025b' ...
           'AlexReis_PCWIN64_MR2025b' ...
           'AlexReis_PCWIN64perso_MR2017a'};
%%% END: TO ADAPT BY USER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
dataFiles   = zol.getFileData(['tex_pdf' filesep 'data' filesep fileDir{1} filesep]);

for i = 1:numel(dataFiles)
    N4_ = []; D4_ = [];
    for ifile = 1:numel(fileDir)
        N4 = []; D4 = []; ZE4 = []; PO4 = []; 
        %
        dataDir     = ['tex_pdf' filesep 'data' filesep fileDir{ifile} filesep];
        dataFiles   = zol.getFileData(dataDir);
        load([dataDir dataFiles{i}])
        for icas = 1:numel(Z4)
            if isempty(NUM{icas})
                [n4,d4] = zol.get_numden(Z4{icas}); 
            else
                n4      = NUM{icas};
                d4      = DEN{icas};
            end
            N4  = [N4 [zeros(robj-numel(n4),1); n4]];
            D4  = [D4 [zeros(robj-numel(d4),1); d4]];
            %ZE4 = [ZE4 [NaN*ones(robj-numel(ZER{j}),1); ZER{j}]];
            %PO4 = [PO4 [NaN*ones(robj-numel(POL{j}),1); POL{j}]];
        end
        N4_(:,:,ifile) = N4;
        D4_(:,:,ifile) = D4;
    end
    %
    degree = robj:-1:0;
    %
    sz  = size(N4_);
    r   = sz(1)-1;
    figure
    for li = 1:sz(1)
        for co = 1:sz(2)
            width  = 1/(sz(2)+1);
            offset = co*width-.5;
            %
            subplot(221), hold on, grid on, axis tight
            zol.boxplot(degree(li)+offset,abs(squeeze(real(N4_(li,co,:)))),width,col(co,:));
            set(gca,'YScale','log')
            ylabel('Real part'), xlabel('Power degree')
            title('Numerator','FontSize',24)
            subplot(223), hold on, grid on, axis tight
            zol.boxplot(degree(li)+offset,abs(squeeze(imag(N4_(li,co,:)))),width,col(co,:));
            set(gca,'YScale','log')
            ylabel('Imaginary part'), xlabel('Power degree')
            %
            subplot(222), hold on, grid on, axis tight
            zol.boxplot(degree(li)+offset,abs(squeeze(real(D4_(li,co,:)))),width,col(co,:));
            set(gca,'YScale','log')
            ylabel('Real part'), xlabel('Power degree')
            title('Denominator','FontSize',24)
            subplot(224), hold on, grid on, axis tight
            zol.boxplot(degree(li)+offset,abs(squeeze(imag(D4_(li,co,:)))),width,col(co,:));
            set(gca,'YScale','log')
            ylabel('Imaginary part'), xlabel('Power degree')
            %
            sgtitle(['Approximation $r=' num2str(r) '$'],'FontSize',28)
        end
    end
    legend(name)
    zol.figSavePDF(['tex_pdf/figures/1a/num_den_r' num2str(r)])
    drawnow
    close all

    % %
    % figure
    % subplot(221), hold on
    % bar(degree,abs(real(N4.')))
    % set(gca,'YScale','log')
    % ylabel('Real part'), xlabel('Power degree')
    % title('Numerator','FontSize',24)
    % subplot(223), hold on
    % bar(degree,abs(imag(N4.')))
    % set(gca,'YScale','log')
    % ylabel('Imaginary part'), xlabel('Power degree')
    % legend(name)
    % %
    % subplot(222), hold on
    % bar(degree,abs(real(D4.')))
    % set(gca,'YScale','log')
    % ylabel('Real part'), xlabel('Power degree')
    % title('Denominator','FontSize',24)
    % subplot(224), hold on
    % bar(degree,abs(imag(D4.')))
    % set(gca,'YScale','log')
    % ylabel('Imaginary part'), xlabel('Power degree')
end
license('inuse')


