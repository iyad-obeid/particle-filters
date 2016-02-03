clc
clear all
close all

dt = 1e-3;
tmax = 30;
t = (0:dt:tmax-dt)';
nNeurons = 1;

f = 0.25;
s = pi/2*sin(2*pi*f*t) + pi/2;


mu_start = pi/3*ones(1,nNeurons);            % mean   0 < mu < pi
mu_end = pi/2*ones(1,nNeurons);
sigma_squared_start = 0.02*ones(1, nNeurons);  % variance
sigma_squared_end   = 0.02*ones(1, nNeurons);
rho_start  = 20*ones(1, nNeurons);              % scaling factor
rho_end = 20*ones(1, nNeurons);
sigma_squared = zeros(length(t), nNeurons);
mu  = zeros(length(t), nNeurons);
rho = zeros(length(t), nNeurons);

for k = 1:nNeurons
    sigma_squared(:, k) = linspace(sigma_squared_start(1, k), sigma_squared_end(1, k), length(t));
    mu(:, k) = linspace(mu_start(1, k), mu_end(1, k), length(t));
    rho(:, k) = linspace(rho_start(1, k), rho_end(1, k), length(t));
end

alpha = ((mu/pi).^2-(mu/pi).^3-(mu/pi).*sigma_squared)./sigma_squared;
beta = ((mu/pi) - 2*(mu/pi).^2 + (mu/pi).^3-sigma_squared)./sigma_squared;
lambda = zeros(length(t), nNeurons);

for k = 1:length(t)
    lambda(k,:) = rho(k,:).*betapdf(s(k,1)/pi, alpha(k,:), beta(k,:));
end

N = poissrnd(lambda*dt);
N = N>0;

clear sigma_squared_end mu_end rho_end c k f
save input_data

figure, plot(t, s(:, 1), t, mu(:, 1), t, s.*N(:, 1), 'ro'), ...
    xlabel('time (s)'), ylabel('\theta_{p}'),...
    xlim([0 tmax]), ylim([0.01 pi]), grid on;
% s = 0:0.01:1;
% figure, plot(s, rho(1,1)*betapdf(s, alpha(1,1), beta(1,1)), ...
%     s, rho(end,1)*betapdf(s, alpha(end,1), beta(end,1)));
