% Loewner algorithm (simple version) 
% Author: C. Poussot-Vassal [MOR Digital Systems & ONERA]
% 
% Syntax
% [hr,info] = zol.loewner(la,mu,W,V)
%  
% Input arguments
%  - la : interpolation points (k x 1, complex)
%  - mu : interpolation points (q x 1, complex)
%  - W  : data evaluated at points "la" (1 x k, complex)
%  - V  : data evaluated at points "mu" (1 x q, complex)
% 
% Output arguments
%  - hr   : approximation model (handle function)
%  - info : structure with informations about the Loewner world
%    * r    : rational order (integer)
%    * la   : lambda's column interpolation points (k x 1, complex)
%    * mu   : mu's row interpolation points (q x 1, complex)
%    * sv   : normalized singular values of [LL SS] (min(q,k), real)
%    * LL   : Loewner matrix (q x k, complex)
%    * SS   : shifted Loewner matrix (q x k, complex)
%    * V,W  : same as input
%    * LA   : Lambda matrix (k x k, complex)
%    * MU   : Mu matrix (q x q, complex)
%    * L,R  : left, right tangential directions, here 1 (data x tangent directions)
%    * lar  : compressed column (right) interpolation points (r x 1, complex)
%    * mur  : compressed row (left) interpolation points (r x 1, complex)
%    * Hr   : compressed Loewner form (state-space, complex)
%               Hr(s)=Cr(sEr-Ar)\Br+Dr, 
%             where (Er,Ar,Br,Cr,Dr) are available in info.Er ...
% 
% Note 
% Sylvester equations 
%   MU*LL-LL*LA = V*R-L*W 
%   MU*SS-SS*LA = MU*V*R-L*W*LA
% may be checked as
%   test1 = info.MU*info.LL - info.LL*info.LA;
%   test2 = info.V*info.R - info.L*info.W;
%   norm(test1-test2) % small
%   test1 = info.MU*info.SS - info.SS*info.LA;
%   test2 = info.MU*info.V*info.R - info.L*info.W*info.LA;
%   norm(test1-test2) % small
% 
% Description
% Loewner rules.
%

function [hr,info] = loewner(la,mu,W,V,opt)

if nargin < 5 || ~isa(opt,'struct')
    D           = 0;
    robj        = inf;
elseif isa(opt,'struct')
    if isfield(opt,'target')
        robj = opt.target;
    else
        robj = inf;
    end
    if isfield(opt,'D')
        D = opt.D;
    else
        D = 0;
    end
end
info = [];
% LL and SS
if ~isempty(intersect(la,mu)) 
    error('Repetition in "la" and "mu"')
end
k   = length(la);
q   = length(mu);
L   = ones(q,1);
R   = ones(1,k);
LL  = zeros(q,k);
SS  = zeros(q,k);
for ii = 1:q
    for jj = 1:k
        num1        = V(ii,:) - W(:,jj);
        num2        = mu(ii)*V(ii,:) - W(:,jj)*la(jj);
        den         = mu(ii)-la(jj);
        LL(ii,jj)   = num1/den;
        SS(ii,jj)   = num2/den;
    end
end
% D-term
if ~norm(D) == 0
    SS  = (SS - L*D*R);
    V   = V - L*D;
    W   = W - D*R;
end    
% Truncate
[L1,S1,~]   = svd([LL,SS],'econ');
[~,~,R2]    = svd([SS',LL']','econ');
sv          = diag(S1)/S1(1,1);
if robj < 1
    r = sum(diag(S1)/S1(1,1)>robj);
else
    r = min(robj,min(k,q));
end
% Projection
Y   = L1(:,1:r);
X   = R2(:,1:r);
Ar  = -Y'*SS*X;
Br  = Y'*V;
Cr  = W*X;
Er  = -Y'*LL*X;
Lr  = Y'*L;
Rr  = R*X;
% Compressed IP
evR = eig(Er\(Ar+Br*Rr));
evL = eig((Ar+Lr*Cr)/(Er));
% Rational function in matrix form
hr  = @(s) Cr*((s*Er-Ar)\Br)+D;
% 
info.r      = length(Ar);
info.la     = la;
info.mu     = mu;
info.sv     = sv;
info.LL     = LL; 
info.SS     = SS;
info.V      = V;
info.W      = W;
info.LA     = diag(info.la);
info.MU     = diag(info.mu);
info.L      = L;
info.R      = R;
info.X      = X;
info.Y      = Y;
info.lar    = evR;
info.mur    = evL;
%
info.Er     = Er;
info.Ar     = Ar;
info.Br     = Br;
info.Cr     = Cr;
info.Dr     = D;
