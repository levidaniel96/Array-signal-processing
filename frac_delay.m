function y = frac_delay( x, delay, fs )

% Input arguments
% x      : discrete input signal
% delay  : signal delay in seconds
% fs     : sampling frequency in Hz
%
% Output arguments
% y      : delayed discrete input signal

tau_disc = delay*fs;   % Delay in samples

% calculation of integer delay (quantizatin)
res = 500;             % precision of the quantization.
temp = round(tau_disc*res);
tau_num = temp/gcd(temp,res);
tau_den = res/gcd(temp,res);

y = sub_delay( x, tau_num, tau_den ); 

end


function y = sub_delay(x,N,D)
%
% Delay of a time series, by sub-sample.
%
% Usage: y = sub_delay(x,N,D)
%
%
% Input: x - time series to be delayed
%        N/D - Delay amount.
%

z = interp(x,D);
zz = delay(z,N);
y = decimate(zz,D,'FIR');
end


function x_delayed = delay(x,D)
%
% Usage:  x_delayed = delay(x,D);
%
% x - the original signal
% D - amount of delay;
% if D > 0 there are D zeros in the begining and the last D samples are lost.
% if D < 0 the first D samples are lost and there are D zeros in the end.
%
% x_delayed - is the delayed version of x;

tr = 0;
[m,n] = size(x);
if m > n
  x = x';
  tr = 1;
end;

x_delayed = zeros(size(x));
	
if D > 0
  
  x_delayed(:,D+1:length(x)) = x(:,1:length(x)-D);

elseif D <= 0
  
  x_delayed(:,1:length(x)-abs(D)) = x(:,abs(D)+1:length(x));

end;

if tr, x_delayed = x_delayed'; end;

end
