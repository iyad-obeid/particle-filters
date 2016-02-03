%clc
clear all
close all
load kinematics

tmax = 30;
t = (0:dt:tmax-dt)';
nNeurons = 30;

% f = 0.2;
% s = pi/2*sin(2*pi*f*t) + pi/2;
s = (kinematics(1,1:length(t)) + abs(min(kinematics(1,:))))';

theta_p_start = pi*rand(1,nNeurons);
flag = true;
while flag
    v = 0.2*randn(1, nNeurons);
    if (theta_p_start+v) > 0 & (theta_p_start+v) < pi
        theta_p_end = theta_p_start+v;
        flag = false;
    end
end

% 
% gamma_start = 5*ones(1, nNeurons);         % amplitude
% gamma_end   = gamma_start;
% beta_start  = gamma_start;         % offset
% beta_end = beta_start;

gamma_start = 20*ones(1, nNeurons);         % amplitude
gamma_end   = gamma_start;
beta_start  = gamma_start;         % offset
beta_end = beta_start;

% gamma_start = 20*rand(1, nNeurons);         % amplitude
% gamma_end   = 20*rand(1, nNeurons);
% beta_start  = gamma_start;         % offset
% beta_end = gamma_end;

gamma = zeros(length(t), nNeurons);
theta_p  = zeros(length(t), nNeurons);
beta = zeros(length(t), nNeurons);

for k = 1:nNeurons
    gamma(:, k) = linspace(gamma_start(1, k), gamma_end(1, k), length(t));
    theta_p(:, k) = linspace(theta_p_start(1, k), theta_p_end(1, k), length(t));
    beta(:, k) = linspace(beta_start(1, k), beta_end(1, k), length(t));
end

lambda = zeros(length(t), nNeurons);
for k = 1:length(t)
    lambda(k,:) = beta(k,:) + gamma(k,:).*cos(2*(s(k,1) - theta_p(k,:)));
end

N = poissrnd(lambda*dt);
N = N>0;

save input_data


% figure,
% plot(t, s(:, 1), t, theta_p(:, 1), t, s.*N(:, 1), 'ro'), ...
%     xlabel('time (s)'), ylabel('\theta_{true}'),...
%     xlim([0 tmax]), ylim([0.1 pi]), grid on

disp(gamma_start);
disp(sum(N));
% 
figure, plot(t,s)
for k = 1:nNeurons
    hold on, plot(t, theta_p(:,k),'k'), grid on, ylim([0 pi]), 
    xlabel('Time (s)'), ylabel('Preferred Directions');
end

% figure,
% for k = 1:nNeurons
%     hold on, plot(t, gamma(:,k),'k'), grid on, ylim([0 20]), 
%     xlabel('Time (s)'), ylabel('Maximum Firing Rates (Hz)');
% end