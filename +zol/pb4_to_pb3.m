% Zolotarev 4th to 3rd problem
% Author: C. Poussot-Vassal [MOR Digital Systems & ONERA]
% 
% Syntax
% [hz,p,sig,tau] = zol.pb4_to_pb3(h,pts,val)
%  
% Input arguments
%  - h   : Z4 rational approxiomation (handle function)
%  - pts : interpolation points (complex vector)
%  - val : data evaluated at points "pts" (complex vector)
% 
% Output arguments
%  - hz  : Z3 rational approxiomation (handle function)
%  - p   : (1-sig)/(1+sig)
%  - sig : approximation value (real)
%  - tau : H-inf norm error of the Z4 rational approximation  
% 
% Description
% Translate Z4 to Z3
%

function [hz,p,sig,tau] = pb4_to_pb3(h,pts,val)

N       = numel(pts);
val_h   = zeros(N,1);
for ii = 1:N
    val_h(ii,1) = h(pts(ii));
end
tau = norm(val_h(:)-val(:),inf);
sig = (tau/(1+sqrt(1-tau^2)))^2;
p   = (1-sig)/(1+sig);
hz  = @(z) sqrt(sig)*(p+h(z))./(p-h(z));
