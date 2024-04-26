
% ASP lab
%
% Exercise: Half-power Beamwidth
%
% Template
%
% H. Loellmann (June 2017)


Ndl_set = [1000:-0.001:1];          % values for N d/lambda
theta_deg =[2.5 5 10 20 30 45 90];  % steering angles in deg
theta_rad = theta_deg*pi/180;       % steering angles in rad

% Half-power bandwidth in the theta-space for endfire
c0 = 0.443;
Theta_H_end = 2*acos( 1 - c0./Ndl_set);  % rad
Theta_H_end_deg = Theta_H_end*180/pi;    % degrees

figure(2)
clf
loglog( Ndl_set, abs(Theta_H_end_deg),'-b')
set(gca,'Fontsize',12)
axis([ Ndl_set(end) Ndl_set(1) 0.04 150])
xlabel('{\itN d} / \lambda','Fontsize',12)
ylabel('3-dB BW in degrees','Fontsize',12)
text(12,35,'Endfire','Fontsize',12)
box on
grid on
hold on
title('halp-power beamwidth for the scan limit')
% Half-power bandwidth (in deg) for scan limit

% insert your code here .....
c0 = 0.443;
Theta_H_end =acos( 1 -  2*c0./Ndl_set);  % rad
Theta_H_end_deg = Theta_H_end*180/pi;    % degrees
loglog( Ndl_set, abs(Theta_H_end_deg),'-r')
text(5,20,'scan limit','Fontsize',12)

% Half-power bandwidth for different angles
bw = zeros(1,length(theta_rad));
figure;
for k = 1 : length(theta_rad)
    %for l = 1 : length(Ndl_set)
        % insert your code here ....
        Theta_H = acos( cos(theta_rad(k)) - c0./Ndl_set)-acos( cos(theta_rad(k)) +c0./Ndl_set);  % rad
        Theta_H_end(k,:) = Theta_H*180/pi;    % degrees
    %end
   
    % insert your code here....
    loglog( Ndl_set, abs(Theta_H_end(k,:)),'-r')
    set(gca,'Fontsize',12)
    axis([ Ndl_set(end) Ndl_set(1) 0.04 150])
    xlabel('{\itN d} / \lambda','Fontsize',12)
    ylabel('3-dB BW in degrees','Fontsize',12)
    box on
    grid on
    hold on
    rightpoint=abs(Theta_H_end(k,1));
    text(1040,rightpoint,[num2str(theta_deg(k)),'^o'],'fontsize',10);
    title('halp-power beamwidth for diffrent steering angles')
end

