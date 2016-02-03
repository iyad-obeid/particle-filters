clc
clear all
close all

dt = 1e-3;
tmax = 10;
t = (0:dt:tmax-dt)';
nNeurons = 25;

f = 0.1;
s = pi/2*sin(2*pi*f*t) + pi/2;

theta_p_start = linspace(0,pi,nNeurons);
theta_p_end = theta_p_start.*ones(1, nNeurons);
gamma_start = 10*ones(1, nNeurons);         % amplitude
gamma_end   = 00*ones(1, nNeurons);
beta_start  = 10*ones(1, nNeurons);         % offset
beta_end = 10*ones(1, nNeurons);
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

% for k = 1:nNeurons
% figure(k), 
% set(gcf, 'WindowStyle', 'docked'),
% plot(t, s(:, 1), t, theta_p(:, k), t, s.*N(:, k), 'ro'), ...
%     xlabel('time (s)'), ylabel('\theta_{true}'),...
%     xlim([0 tmax]), ylim([0.1 2*pi]), grid on, hold on, ...
%     if theta_p(:, k) < pi/2
%     plot([0 tmax],[(theta_p(:, k)+pi/2) (theta_p(:, k)+pi/2)]);
%     end
%     if theta_p(:, k) > pi/2
%         plot([0 tmax],[(theta_p(:, k)-pi/2) (theta_p(:, k)-pi/2)]);
%     end
% end
    


%     figure, plot(t, s(:, 1), t, theta_p(:, 1), t, s.*N(:, 1), 'ro'), ...
%     xlabel('time (s)'), ylabel('\theta_{true}'),...
%     xlim([0 tmax]), ylim([0.1 2*pi]), grid on, hold on, ...
%     if theta_p(:, 1) < pi/2
%     plot([0 tmax],[(theta_p(:, 1)+pi/2) (theta_p(:, 1)+pi/2)]);
%     end
%     if theta_p(:, 1) > pi/2
%         plot([0 tmax],[(theta_p(:, 1)-pi/2) (theta_p(:, 1)-pi/2)]);
%     end
    
    
%     for k = 2:nNeurons
%         hold on, plot(t, theta_p(:,k))
%     end

    % s_L = 0:0.01:pi;
    % L = beta(1,1) + gamma(1,1).*cos(2*(s_L - theta_p(1,1)));
    % figure, plot(s_L, L), xlim([0 pi]), ...
    %     xlabel('\theta'), ylabel('firing rate (Hz)'), grid on;
