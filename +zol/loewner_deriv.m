
function [hr,info] = loewner_deriv(val,pts,ROBJ,D,shift)
    
    info    = [];
    shift   = repmat([shift; -shift],length(val)/2,1);
    la      = pts;
    mu      = la.';
    W(1,:)  = val+shift;
    V       = W.';
    
    %%% LL and SS
    k   = length(la);
    q   = length(mu);
    L   = ones(q,1);
    R   = ones(1,k);
    LL  = zeros(q,k);
    SS  = zeros(q,k);
    for ii = 1:q
        for jj = 1:k
            if ii==jj
                LL(ii,jj)   = 0;
                SS(ii,jj)   = W(1,ii);
            else 
                num1        = V(ii,:) - W(:,jj);
                num2        = mu(ii)*V(ii,:) - W(:,jj)*la(jj);
                den         = mu(ii)-la(jj);
                LL(ii,jj)   = num1/den;
                SS(ii,jj)   = num2/den;
            end
        end
    end
    % D-term
    % if ~norm(D) == 0
    %     SS  = (SS - L*D*R);
    %     V   = V - L*D;
    %     W   = W - D*R;
    % end

    % %%% Go real
    % TOL_CC = 1e-13;
    % isCC = false;
    % if (abs(sum(imag(la_.')))<TOL_CC) && ...
    %    (abs(sum(imag(mu_.')))<TOL_CC) && ...
    %    (q==k)
    %     isCC = true
    %     warning('Go real')
    % end
    % if isCC
    %     J0  = (1/sqrt(2))*[1 1i; 1 -1i];
    %     J   = [];
    %     kk  = 1;
    %     while length(J) < length(la_)
    %         if imag(la_(kk)) == 0
    %             J   = blkdiag(J,1);
    %             kk  = kk + 1;
    %         else
    %             J   = blkdiag(J,J0);
    %             kk  = kk + 2;
    %         end
    %     end
    %     LL  = real(J'*LL*J);
    %     SS  = real(J'*SS*J);
    %     V   = real(J'*V);
    %     W   = real(W*J);
    % end
    
    %%% Truncate
    [L1,S1,~]   = svd([LL,SS],'econ');
    [~,~,R2]    = svd([SS',LL']','econ');
    sv          = diag(S1)/S1(1,1);
    if ROBJ < 1
        r = sum(diag(S1)/S1(1,1)>ROBJ);
    else
        r = min(ROBJ,q);
    end
    %%% Projection
    Y   = L1(:,1:r);
    X   = R2(:,1:r);
    Ar  = -Y'*SS*X;
    Br  = Y'*V;
    Cr  = W*X;
    Er  = -Y'*LL*X;
    Lr  = Y'*L;
    Rr  = R*X;
    %%% Compressed IP
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
end
