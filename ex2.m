%% ex2 -  Non-uniform Weighting
% parameters
Sampels=1000;
Vec_N=11;
N=11;
w_n=1/N;
n=0:N-1;
lambda=1;
d=lambda/2;
n0=(N-1)/2;
n_m_n0=n-n0;
%% u space 
w1=w_n.*ones(N,1);
w2=sin(pi/(2*N)).*cos(pi*n_m_n0/N).';
u=linspace(0,1,Sampels);
psi=2*pi*d/lambda*u;
B_1=ULA(psi,w1,n_m_n0);
B_2=ULA(psi,w2,n_m_n0);
plot(u,10*log10(abs(B_1.^2)))
hold on
plot(u,10*log10(abs(B_2.^2)))
legend('uniform Weighting',' Non-uniform Weighting')
ylim([-100,0])
ylabel('powerpattern[dB]')
xlabel('u')
title('Non-uniform Weighting comprasion')