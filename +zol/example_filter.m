% List of Zolotarev examples 
% Author: C. Poussot-Vassal [MOR Digital Systems & ONERA]
% 
% Syntax
% [pts,val,data] = zol.example_filter(name)
%  
% Input arguments
%  - name : 
% 
% Output arguments
%  - pts  : interpolation points (complex vector)
%  - val  : data evaluated at points "pts" (+/- 1 vector)
%  - data : structure of the topology
% 
% Description
% List of exmaples where E = -1 (left) and F = +1 (right)
% 

function [pts,val,info] = example_filter(wn,amp)

Xlim    = [-1 1]*1.5;
Ylim    = [-1 1]*1.5;
E       = [];
F       = [];
maxAmp  = max(amp);
minAmp  = min(amp);
k       = 1;
while k < length(wn)
    if amp(k) == minAmp
        E = [E exp(1i*linspace(wn(k),wn(k+1),50)*pi)];
    elseif amp(k) == maxAmp
        F = [F exp(1i*linspace(wn(k),wn(k+1),100)*pi)];
    end
    k = k + 2;
end
E   = [E conj(E)];
F   = [F conj(F)];

%%% Data
la          = E(:); 
n1          = length(la);
mu          = F(:); 
n2          = length(mu);
W(1:n1,1)   = minAmp;
V(1:n2,1)   = maxAmp;
pts         = [la; mu];
val         = [W; V];

% 
info.E      = E;
info.F      = F;
info.bnd    = [minAmp maxAmp]; % [E F]
info.Xlim   = Xlim;
info.Ylim   = Ylim;


% %%
% function x = zol.chebspace(a,b,n)
% 
%     x = zeros(1,n);
%     for k = 0:n-1
%         x(1,k+1) = cos((k+.5)*pi/n);
%     end
% 
%     x = .5*(a+b)+.5*(b-a)*x;
% 
% end