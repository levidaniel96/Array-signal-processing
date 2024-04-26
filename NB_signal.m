
function [y,h] = NB_signal( fs, f0, N)

% Input arguments:
%
% fs : Sampling frequency in Hz
% f0 : Carrier frequency in Hz
% N  : Signal length in samples
%
% Output arguements:
% y  : real-valued NB signals (bandlimited white noise)
% Author: S. Gannot, BIU

if nargin == 0
    fs = 16e3;
    f0 = fs/3;
    N = 2^10;    
end

% Lowpass filter
L = 151; % FIR order
dw = pi/500; % bandwidth
h = fir1(L,dw/pi);
j=sqrt(-1);

w0_d = 2*pi*f0/fs; % discrete center frequency

carrier = exp(j*w0_d*[0:N-1]);
s = randn(1,N);
s_nb = filter(h,1,s); % narrow band white noise
%sa = s_nb+j*hilbert(s_nb); % analytical signal no carrier
sa = hilbert(s_nb); % analytical signal no carrier
sa_c = sa.*carrier; % analytical signal with carrier

y = real(sa_c); 