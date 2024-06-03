%% ex1 - ULA 
% basic parameters
Sampels=1000;
Vec_N=[10,11];
for idx_N=1:length(Vec_N)
    % ULA parameters
    N=Vec_N(idx_N); % num of sensors
    w_n=1/N; % ULA weights
    w=w_n.*ones(N,1);
    n=0:N-1;
    lambda=1; 
    d=lambda/2;
    n0=ceil((N-1)/2);
    n_m_n0=n-n0;
%% a)

    %% Kz space 
    figure;
    sgtitle(['ULA Beampattern for N=',num2str(N)])
    subplot(4,1,1)
    Kz=linspace(-1,1,Sampels)*2*pi/lambda;
    psi=-Kz*d;
    B_Kz=ULA(psi,w,n_m_n0);
    plot(Kz,abs(B_Kz))
    title('Kz space')
        
    %% psi space 
    subplot(4,1,2)
    psi=linspace(-1,1,Sampels)*2*pi*d/lambda;
    B_psi=ULA(psi,w,n_m_n0);
    plot(psi,abs(B_psi))
    title('\psi space')

    %% u space 
    subplot(4,1,3)
    u=linspace(-1,1,Sampels);
    psi=2*pi*d/lambda*u;
    B_u=ULA(psi,w,n_m_n0);
    plot(u,abs(B_u))
    title('u space')
    %% theta space 
    subplot(4,1,4)
    theta=linspace(0,pi,Sampels);
    psi=2*pi*d/lambda*cos(theta);
    B_theta=ULA(psi,w,n_m_n0);
    plot(theta,abs(B_theta))
    title('\theta space')
    xlim([0 pi])
    
%% b)
    figure
    theta=linspace(0,pi,Sampels);
    psi=2*pi*d/lambda*cos(theta);
    B_theta=ULA(psi,w,n_m_n0);
    polardb(theta,10*log10(abs(B_theta.^2)),-40)
    title(['Powerpattern for N =',num2str(N)])
end
%% c)
% parameters 
N=10;
n=0:N-1;
n0=ceil((N-1)/2);
n_m_n0=n-n0;
% broadside
w_n=1/N;
w=w_n.*ones(N,1);
theta=linspace(0,pi,Sampels);
d_lambda=linspace(0.001,1,Sampels)';
for idx=1:length(d_lambda)
    psi=2*pi*d_lambda(idx)*cos(theta);
    B(idx,:)=ULA(psi,w,n_m_n0);     
end
figure
imagesc(theta,d_lambda,10*log10(abs(B.^2)))
caxis([-40,0])
ylabel('d/\lambda')
xlabel('\Theta')
title(['powerpattern in the \Theta space broadside'])
% endfire
theta=linspace(0,pi,Sampels);
d_lambda=linspace(0.001,1,Sampels)';
ang=0;
for idx=1:length(d_lambda)
    psi=2*pi*d_lambda(idx)*cos(theta);
    psi2=2*pi*d_lambda(idx)*cos(ang);
    w=1/N*exp(1i*(n_m_n0)'*psi2);
    B(idx,:)=ULA(psi,w,n_m_n0);     
end
figure
imagesc(theta,d_lambda,10*log10(abs(B.^2)))
caxis([-40,0])
ylabel('d/\lambda')
xlabel('\Theta')
title(['powerpattern in the \Theta space endfire'])

%% d)
% parameters 
N=10;
n=0:N-1;
n0=ceil((N-1)/2);
n_m_n0=n-n0;
w_n=1/N;
theta=linspace(0,pi,Sampels);
d_lambda=linspace(0.001,1,Sampels)';
ang_vec=[0,30,60,90];

for ang_idx=1:length(ang_vec)
    ang=deg2rad(ang_vec(ang_idx));
    for idx=1:length(d_lambda)
        w=w_n.*exp(1i*n_m_n0*2*pi*d_lambda(idx)*cos(ang));
        psi=2*pi*d_lambda(idx)*cos(theta);
        B(idx,:)=ULA(psi,w.',n_m_n0);     
    end
    
subplot(2,2,ang_idx)
imagesc(rad2deg(theta),d_lambda,10*log10(abs(B.^2)))
caxis([-40,0])
ylabel('d/\lambda')
xlabel('\Theta')
title(['angle= ',int2str(ang_vec(ang_idx))])
end

%% e)
figure
N=11;
w_n=1/N;
n=0:N-1;
lambda=1;
d=lambda/2;
n0=ceil((N-1)/2);
n_m_n0=n-n0;
w=w_n.*ones(N,1);
figure;
psi=linspace(-1,1,Sampels)*2*pi*d/lambda;
% all sensor works
B_psi=ULA(psi,w,n_m_n0);
plot(psi,abs(B_psi))
hold on 
% failure sensor case
w([3,5,6])=0;
B_psi=ULA(psi,w,n_m_n0);
plot(psi,abs(B_psi.^2))
legend('all sensors works','sensors 3,5,6 doesnt works')