%% ex4
% parameters
Sampels=1000;
N=10;
n=0:N-1;
n0=(N-1)/2;
n_m_n0=n-n0;
sigma_w_2=1;
SNR=70;
sigma_1_2=sigma_w_2*10^(SNR/10);
v_s=1/N.*ones(N,1);
d_lambda=1/2;
V_u=[0.3,0.004];
for i=1:length(V_u)
    u1=V_u(i);
    psi=2*pi*d_lambda*u1;
    v_1=exp(1i*n_m_n0*psi).';
    Sn=sigma_w_2*eye(N)+sigma_1_2*v_1*v_1';

    w_MVDR=v_s'*inv(Sn)/(v_s'*inv(Sn)*v_s);

    %% u space 
    u=linspace(-1,1,Sampels);
    psi=2*pi*d_lambda*u;
    B_1=ULA(psi,w_MVDR',n_m_n0);
    plot(u,10*log10(abs(B_1.^2)))
    hold on 
end
ylabel('powerpattern[dB]')
xlabel('u')
legend('u1=0.3','u1=0.004')

%% b
Sampels=1000;
N=10;
n=0:N-1;
n0=(N-1)/2;
n_m_n0=n-n0;
sigma_w_2=1;
SNR_vec=[70,0];
for SNR_idx=1:length(SNR_vec)
    SNR=SNR_vec(SNR_idx);
    sigma_1_2=sigma_w_2*10^(SNR/10);
    v_s=1/N.*ones(N,1);
    d_lambda=1/2;
    u1=linspace(0.001,0.5,Sampels);
    u=linspace(-1,1,Sampels);

    for idx=1:length(u1)
        psi=2*pi*d_lambda*u1(idx);
        v_1=exp(1i*n_m_n0*psi).';
        Sn=sigma_w_2*eye(N)+sigma_1_2*v_1*v_1';
        w_MVDR=v_s'*inv(Sn)/(v_s'*inv(Sn)*v_s);
        psi=2*pi*d_lambda*u;
        B(idx,:)=ULA(psi,w_MVDR',n_m_n0);     
    end
    figure
    imagesc(u,u1,10*log10(abs(B.^2)))
    caxis([-100,45])
    ylabel('u_1')
    xlabel('u')
    colorbar
    title(['powerpattern for varying interferer directions for SNR=',num2str(SNR)])
end