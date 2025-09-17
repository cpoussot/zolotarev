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
