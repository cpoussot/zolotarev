% Main interface for chirp signal generation
% Author: C. Poussot-Vassal [MOR Digital Systems / Onera]
%
% Description
% Computes a chirp signal sweeping over FBND frequency range.
%
% Syntax
%  [uc,tc,info] = insapack.chirp(Ns,Ts,FBND,REV,TYPE,SHOW)
% 
% Input arguments
%  - Ns   : number of samples (integer)
%  - Ts   : sampling period [s] (>0 real)
%  - FBND : frequency range [Hz] (2x1 real vector)
% 
% Output arguments
%  - uc   : chirp signal in [-1,1]
%  - tc   : time samples
% 

function [uc,tc] = chirp(Ns,Ts,FBND)

% Nyquist check
Fs  = 1/Ts;
if max(FBND) >= Fs/2
    error('Nyquist frequency issue: check "max(fBnd) < 1/(2Ts)".')
end
%
tc  = (0:Ns-1)*Ts;
f0  = min(FBND);
f1  = max(FBND);
t1  = tc(end);
% Set p=1 for linear, 2 for quadratic, 3 for logarithmic
p       = 1;
beta    = (f1-f0).*(t1.^(-p));
uc      = sin(2*pi * ( beta./(1+p).*(tc.^(1+p)) + f0.*tc));
% % FFT
% L       = numel(uc);
% FTuc    = fft(uc)/L;
% fc      = linspace(0,1,L)*Fs;
