%% ex5
% Basic parameters
M=10;
d=0.05;
c=342;
Samples=1024;
fs=8000;
f0=7/3*10^3;
lambda = c/f0;
% Create a NB signal
[x,h]=NB_signal(fs, f0, Samples);

ang_vec = 0:5:90;
for ang_idx = 1:length(ang_vec)
    ang=ang_vec(ang_idx);
    delay=-d/c*cos(deg2rad(ang));
    for m=0:M-1
        y(m+1,:) = frac_delay( x, m*delay, fs );
        x_tilda(m+1,:) = hilbert(y(m+1,:));

    end        
    % Beamforming in base-band
    theta_deg_vec = [0,90];
    for theta_idx=1:length(theta_deg_vec)
        theta_deg=theta_deg_vec(theta_idx);
        w= exp( 1j*(0:M-1)'*2*pi*d/lambda*cos(deg2rad(theta_deg)) )/M;
        Y  = w'*x_tilda;
        OIR(theta_idx,ang_idx)  = var(Y)/var(x);
    end
end

% ULA parameters
w_n=1/M;
w=w_n.*ones(M,1);
n=0:M-1;
n0=ceil((M-1)/2);
n_m_n0=n-n0;
%% theta space
% endfire
w_n=1/M;
w=w_n.*ones(M,1);
theta= 0:5:90;
ang=90;
psi=2*pi*d/lambda*cos(deg2rad(theta));
psi2=2*pi*d/lambda*cos(deg2rad(ang));
w=1/M*exp(1i*(n_m_n0)'*psi2);
B=ULA(psi,w,n_m_n0);
figure(1)
plot(ang_vec,OIR(2,:)/max(OIR(2,:)),'-o')
hold on
grid
xlabel('DOA [degree]')
plot(theta,abs(B.^2))
legend('OIR','beampattern')

% broadside
ang=0;
psi=2*pi*d/lambda*cos(deg2rad(theta));
psi2=2*pi*d/lambda*cos(deg2rad(ang));
w=1/M*exp(1i*(n_m_n0)'*psi2);
B=ULA(psi,w,n_m_n0);

figure(2)
plot(ang_vec,OIR(1,:)/max(OIR(1,:)),'-o')
hold on
grid
xlabel('DOA [degree]')
plot(theta,abs(B.^2))
legend('OIR','beampattern')
