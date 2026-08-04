clearvars; close all; clc; format shorte
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN: TO ADAPT BY USER
% Add Zolotarev Loewner package
addpath('/Users/charles/Documents/GIT/zolotarev')
% Add AAA package
addpath('/Users/charles/Documents/GIT/_others/chebfun')
% Your name: e.g. for Charles Poussot-Vassal: 'cpv'
nameUsr = 'cpv';
%%% END: TO ADAPT BY USER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Create folder
fileDir     = [nameUsr '_' computer];
fileData    = ['tex_pdf' filesep 'data' filesep fileDir filesep];
if ~exist(fileData); mkdir(fileData); end

%%% Choose case and paramters for AAA
AAAparam    = {... % ",'sign',1,'lawson',0"; ... 
               ",'sign',1,'damping',.95,'lawson',200"; ...
               ",'sign',1,'damping',.95,'lawson',1000"};

%%% Style & slack variables
AAAname     = {... %'AAA-sign-L0'; ...
               'AAA-sign-L200'; ...
               'AAA-sign-L1000'};
mw          = 20;
lw          = 4;
col         = parula(10);
col(1,:)    = [0 0 0];
col(2,:)    = [1 0 0];
mark        = {'s' 'o' '<' '>' '^' 'v'};

%%% Define Zolotarev topology
[pts,val,data]  = zol.example('1a');
for robj = 3:10
    k       = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Optimal solution
    k       = k+1;
    Z4{k}   = data.z4_opt{robj};
    SIG{k}  = data.sig_opt{robj};
    NUM{k}  = data.num_opt{robj};
    DEN{k}  = data.den_opt{robj};
    ZER{k}  = data.zer_opt{robj};
    POL{k}  = data.pol_opt{robj};
    name{k} = ['Opt., $\sigma_{' num2str(robj) '}=' num2str(SIG{k},2) '$'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Loewner
    k                   = k+1;
    [la,mu,W,V]         = zol.example2data(pts,val,data);
    opt.target          = robj;
    [h4_loe,info]       = zol.loewner(la,mu,W,V,opt);
    % [n4_loe,d4_loe,M]   = zol.get_numden(h4_loe); 
    [h3_loe,~,sig_loe]  = zol.pb4_to_pb3(h4_loe,pts,val);
    H_loe               = dss(info.Ar,info.Br,info.Cr,0,info.Er,0);
    Z4{k}               = h4_loe;
    SIG{k}              = sig_loe;
    NUM{k}              = [];%n4_loe
    DEN{k}              = [];%d4_loe
    ZER{k}              = sort(zero(H_loe)); %eig([info.Ar info.Br;info.Cr 0],blkdiag(info.Er,0));
    POL{k}              = sort(eig(H_loe));  %eig(info.Ar,info.Er);
    name{k}             = ['LF, $\sigma_{' num2str(robj)  '}=' num2str(SIG{k},2) '$'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% AAA
    for ii = 1:length(AAAparam)
        k                   = k+1;
        AAAparam_i          = AAAparam{ii};
        eval(['[h4_aaa,h4_aaa_pol,h4_aaa_res,h4_aaa_zer,zj,fj,wj,errvec] = aaa(val,pts,"degree",robj' AAAparam_i{1} ');'])
        % [n4_aaa,d4_aaa,M]   = zol.get_numden(h4_aaa); 
        [h3_aaa,~,sig_aaa]  = zol.pb4_to_pb3(h4_aaa,pts,val);
        %
        Z4{k}               = h4_aaa;
        SIG{k}              = sig_aaa;
        NUM{k}              = [];%n4_aaa;
        DEN{k}              = [];%d4_aaa;
        ZER{k}              = sort(h4_aaa_zer);
        POL{k}              = sort(h4_aaa_pol);
        name{k}             = [AAAname{ii} ', $\sigma_{' num2str(robj) '} =' num2str(SIG{k},2) '$'];
    end
    %%% Save all
    save([fileData 'data_r' num2str(robj)], ...
          'AAAparam','AAAname','name','data', ...
          'robj','Z4','SIG','ZER','POL','NUM','DEN', ...
          'zj','fj','wj','errvec')
end
license('inuse')
