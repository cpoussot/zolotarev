clearvars; close all; clc; format shorte
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN: TO ADAPT BY USER
% Add Zolotarev Loewner package
addpath('/Users/charles/Documents/GIT/zolotarev')
% Add AAA package
addpath('/Users/charles/Documents/GIT/_others/chebfun')
% Directory name list
fileDir = {'cpv_MACA64' 'cpv_winos'};
%%% END: TO ADAPT BY USER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
%dataFiles   = zol.getFileData(['tex_pdf' filesep 'data' filesep fileDir{1} filesep]);
%syms z
for i = 2%:3%1:numel(dataFiles)
    for ii = 1:numel(fileDir)
        N4 = []; D4 = []; ZE4 = []; PO4 = []; 
        dataDir     = ['tex_pdf' filesep 'data' filesep fileDir{ii} filesep];
        dataFiles   = zol.getFileData(dataDir);
        load([dataDir dataFiles{i}])
        for j = 1:numel(Z4)
            if isempty(NUM{j})
                [n4,d4] = zol.get_numden(Z4{j}); 
            else
                n4      = NUM{j};
                d4      = DEN{j};
            end
            N4  = [N4 [zeros(robj-numel(n4),1); n4]];
            D4  = [D4 [zeros(robj-numel(d4),1); d4]];
            %ZE4 = [ZE4 [NaN*ones(robj-numel(ZER{j}),1); ZER{j}]];
            %PO4 = [PO4 [NaN*ones(robj-numel(POL{j}),1); POL{j}]];
        end
        N4_(:,:,i) = N4;
        D4_(:,:,i) = D4;
    end
    %
    degree = robj:-1:0;
    %
    figure
    for 
    
    zol.boxplot(3,y,.1,'r')


    %
    figure
    subplot(221), hold on
    bar(degree,abs(real(N4.')))
    set(gca,'YScale','log')
    ylabel('Real part'), xlabel('Power degree')
    title('Numerator','FontSize',24)
    subplot(223), hold on
    bar(degree,abs(imag(N4.')))
    set(gca,'YScale','log')
    ylabel('Imaginary part'), xlabel('Power degree')
    legend(name)
    %
    subplot(222), hold on
    bar(degree,abs(real(D4.')))
    set(gca,'YScale','log')
    ylabel('Real part'), xlabel('Power degree')
    title('Denominator','FontSize',24)
    subplot(224), hold on
    bar(degree,abs(imag(D4.')))
    set(gca,'YScale','log')
    ylabel('Imaginary part'), xlabel('Power degree')
    
end
license('inuse')


