% Convert example data to LF
% Author: C. Poussot-Vassal [MOR Digital Systems & ONERA]
% 
% Syntax
% [la,mu,W,V] = zol.example2data(pts,val,data,interlace)
%  
% Input arguments
%  - pts       : interpolation points (complex vector)
%  - val       : data evaluated at points "pts" (+/- 1 vector)
%  - data      : structure of the topology (see zol.example)
%  - interlace : interlace data (boolean, default true)
%  - CC        : complex conjugate IP (boolean, default false)
% 
% Output arguments
%  - la : interpolation points (k x 1, complex)
%  - mu : interpolation points (q x 1, complex)
%  - W  : data evaluated at points "la" (1 x k, complex)
%  - V  : data evaluated at points "mu" (1 x q, complex)
% 
% Description
% Convert example data to LF
%

function [la,mu,W,V] = example2data(pts,val,data,interlace,CC)

if nargin < 4
    interlace   = true;
    CC          = false;
end
if nargin < 5
    CC  = false;
end

%%% Keep as many points in E and F 
% /!\ this is only to the Zolotarev problem, not a limitation of LF
n1      = length(data.E);
n2      = length(data.F);
ratio   = n1/n2;
if ratio == 1
    Nla = 1;
    Nmu = 1;
elseif ratio < 1
    Nla = 1;
    Nmu = round(1/ratio);
elseif ratio > 1
    Nla = round(ratio);
    Nmu = 1;
end
%
la          = pts(1:n1);
mu          = pts(end-n2+1:end);
W(1,1:n1)   = val(1:n1);
V(1:n2,1)   = val(end-n2+1:end);
la          = la(1:Nla:n1);
mu          = mu(1:Nmu:n2);
W           = W(1,1:Nla:n1);
V           = V(1:Nmu:n2,1);
nn          = floor(min(length(la),length(mu))/2)*2;

%%% Interlace IP
% /!\ this is better for the Zolotarev problem, not a limitation of LF
if interlace 
    for i = 1:2:nn
        la_(i)    = la(i);
        la_(i+1)  = mu(i);
        mu_(i)    = la(i+1);
        mu_(i+1)  = mu(i+1);
        W_(1,i)   = W(1,i);
        W_(1,i+1) = V(i,1);
        V_(i,1)   = W(1,i+1);
        V_(i+1,1) = V(i+1,1);
    end
    la  = la_;
    mu  = mu_;
    V   = V_;
    W   = W_;
end

%%% Complex conjugate 
% /!\ if you wanna have real realization in LF
if CC
    ila = find(imag(la)>0);
    imu = find(imag(mu)>0);
    la  = la(ila);
    mu  = mu(imu);
    W   = W(ila);
    V   = V(imu);
    k   = length(la);
    q   = length(mu);
    %R   = ones(1,2*k);
    %L   = ones(2*q,1);
    kk  = 1;
    for ii = 1:k
        la_(kk)     = la(ii); 
        W_(1,kk)    = W(ii); kk = kk +1;
        la_(kk)     = conj(la(ii));
        W_(1,kk)    = W(ii); kk = kk +1;
    end
    kk  = 1;
    for ii = 1:q
        mu_(kk)     = mu(ii); 
        V_(kk,1)    = V(ii); kk = kk +1;
        mu_(kk)     = conj(mu(ii));
        V_(kk,1)    = V(ii); kk = kk +1;
    end
    la  = la_;
    mu  = mu_;
    V   = V_;
    W   = W_;
end