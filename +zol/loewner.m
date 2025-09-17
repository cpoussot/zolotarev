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
