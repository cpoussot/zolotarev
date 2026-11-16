% List of Zolotarev examples 
% Author: C. Poussot-Vassal [MOR Digital Systems & ONERA]
% 
% Syntax
% [pts,val,data] = zol.example_filter(wn,amp,Nw)
%  
% Input arguments
%  - wn : normalized frequency band
%  - amp: amplitude values min/max
%  - Nw : number of points within each band
% 
% Output arguments
%  - pts  : interpolation points (complex vector)
%  - val  : data evaluated at points "pts" (+/- 1 vector)
%  - data : structure of the topology
% 
% Description
% List of exmaples where E = -1 (left) and F = +1 (right)
% 

function [pts,val,info] = example_filter_smooth(wn,amp,Nw)

    Nw      = 2*floor(Nw/2);
    Xlim    = [-1 1]*1.5;
    Ylim    = [-1 1]*1.5;
    maxAmp  = max(amp);
    minAmp  = min(amp);
    maxNum  = sum(amp==maxAmp);
    minNum  = sum(amp==minAmp);
    Nmult   = (maxNum/minNum);
    k       = 1;
    kk      = 0;
    ll      = 0;
    while k < length(wn)
        if amp(k) == minAmp
            kk          = kk+1;
            E{kk}       = exp(1i*linspace(wn(k),wn(k+1),floor(Nmult*Nw))*pi);
            valE{kk}    = minAmp*ones(1,floor(Nmult*Nw));
        elseif amp(k) == maxAmp
            ll          = ll+1;
            F{ll}       = exp(1i*linspace(wn(k),wn(k+1),Nw)*pi);
            valF{ll}    = maxAmp*ones(1,Nw);
        end
        k = k + 2;
    end
    %%% Complete data
    kk      = 0;
    Ncomp   = 2*floor(Nw/2);
    for k = 1:length(wn)-1
        if amp(k) ~= amp(k+1)
            kk = kk + 1;
            EFcomp{kk}      = exp(1i*chebspace(1.01*wn(k),wn(k+1)*.99,Ncomp)*pi);
            valEFcomp{kk}   = amp(k)*chebspace(1,0,Ncomp)+amp(k+1)*chebspace(0,1,Ncomp);
        end
    end
    
    Ev = []; valEv = [];
    Fv = []; valFv = [];
    for ii = 1:numel(E)
        Ev      = [Ev    E{ii}    EFcomp{ii}(2:2:end)];
        valEv   = [valEv valE{ii} valEFcomp{ii}(2:2:end)];
    end
    while ii<numel(EFcomp)
        Ev = [Ev EFcomp{ii+1}(2:2:end)];
        valEv = [valEv valEFcomp{ii+1}(2:2:end)];
        ii = ii + 1;
    end
    for ii = 1:numel(F)
        Fv      = [Fv    F{ii}    EFcomp{ii}(1:2:end)];
        valFv   = [valFv valF{ii} valEFcomp{ii}(1:2:end)];
    end
    while ii<numel(EFcomp)
        Fv = [Fv EFcomp{ii+1}(1:2:end)];
        valFv = [valFv valEFcomp{ii+1}(1:2:end)];
        ii = ii + 1;
    end
    
    %
    E       = [Ev conj(Ev)];
    F       = [Fv conj(Fv)];
    valE    = [valEv valEv];
    valF    = [valFv valFv];
    
    %%% Data
    la      = E(:);
    mu      = F(:);
    W(:,1)  = valE;%minAmp;
    V(:,1)  = valF;%maxAmp;
    pts     = [la; mu];
    val     = [W; V];
    
    % 
    info.E      = E;
    info.F      = F;
    info.bnd    = [minAmp maxAmp]; % [E F]
    info.Xlim   = Xlim;
    info.Ylim   = Ylim;
end
%%
function x = chebspace(a,b,n)
    x = zeros(1,n);
    for k = 0:n-1
        x(1,k+1) = cos((k+.5)*pi/n);
    end

    x = .5*(a+b)+.5*(b-a)*x;
end