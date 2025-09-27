function [hnum,hden,M] = get_numden(h,SHOW,OPT)

if nargin < 2
    SHOW    = false;
    OPT     = [];
end
if nargin < 3
    OPT     = [];
end

hsym    = sym(h);
[n,d]   = numden(hsym);
tmp_n   = sym2poly(n);
tmp_d   = sym2poly(d);
hnum    = tmp_n.'/tmp_d(1);
hden    = tmp_d.'/tmp_d(1);
syms s
for i = 1:numel(hden)
    M(i,1) = s^(numel(hden)-i);
end

%%% Reshape
k   = max(length(hnum),length(hden));
if length(hnum) < k
    hnum = [zeros(k-length(hnum),1); hnum(:)];
end
if length(hden) < k
    hden = [zeros(k-length(hden),1); hden(:)];
end

%%%
if SHOW
    lw  = 4;
    %vpa(hnum,3)
    %vpa(hden,3)
    figure,
    subplot(211), hold on, grid on, axis tight
    plot(0:k-1,flipud(real(hnum)),'o','MarkerSize',20,'linewidth',lw)
    plot(0:k-1,flipud(real(hden)),'s','MarkerSize',20,'linewidth',lw)

    legend({'Numerator' 'Denominator'})
    if ~isempty(OPT)
        plot(0:k-1,flipud(real(OPT{1})),'kx','MarkerSize',20,'linewidth',lw)
        plot(0:k-1,flipud(real(OPT{2})),'.','MarkerSize',30,'linewidth',lw)
        legend({'Numerator' 'Denominator' 'Numerator (opt.)' 'Denominator (opt.)'},'Location','Best')
    end
    xlabel('$s^k$')
    title('Polynomial coefficients real part')
    subplot(212), hold on, grid on, axis tight
    plot(0:k-1,flipud(imag(hnum)),'o','MarkerSize',20,'linewidth',lw)
    plot(0:k-1,flipud(imag(hden)),'s','MarkerSize',20,'linewidth',lw)
    legend({'Numerator' 'Denominator'})
    if ~isempty(OPT)
        plot(0:k-1,flipud(imag(OPT{1})),'kx','MarkerSize',20,'linewidth',lw)
        plot(0:k-1,flipud(imag(OPT{2})),'.','MarkerSize',30,'linewidth',lw)
        legend({'Numerator' 'Denominator' 'Numerator (opt.)' 'Denominator (opt.)'},'Location','Best')
    end
    title('Polynomial coefficients imaginary part')
end