clc
clear
load kinematics

training_period = 300;
angle = kinematics(1:(training_period/dt), 2);
clear kinematics

N = false(length(angle), 158);
for n = 1:158
    eval(['load neuron', num2str(n)]);
    eval(['temp = N', num2str(n),'(1:length(angle));']);
    N(:, n) = temp';
    eval(['clear N', num2str(n), ' N', num2str(n)]);
end

F_elbow = zeros(size(N));
for n = 1:158
    F_elbow(:,n) = N(:,n).*angle;
end

db = 0.05;
bin_centers = db:db:pi;
angle_hist = hist(angle, bin_centers);

neuron_params(1:158)= struct('pos',0);
pd = zeros(1, 158);
gamma_const = zeros(1, 158);
beta_const = zeros(1, 158);

for n = 1:158
    neuron_params(n).pos = find(F_elbow(:,n));
    neuron_params(n).h = hist(angle(neuron_params(n).pos), bin_centers);
    neuron_params(n).h_norm = (neuron_params(n).h)./angle_hist;
    
    [~, index] = max(neuron_params(n).h_norm);
    pd(n) = bin_centers(index);
    neuron_params(n).max_rate = (neuron_params(n).h(index))/(angle_hist(index)*dt);
    [~, index] = min(neuron_params(n).h_norm);
    neuron_params(n).min_rate = (neuron_params(n).h(index))/((angle_hist(index)*dt));
    
    gamma_const(1, n) = (neuron_params(n).max_rate - neuron_params(n).min_rate)/2;
    beta_const(1, n) = gamma_const(1, n) + neuron_params(n).min_rate;
%     
%     figure(n), set(gcf,'WindowStyle','docked'), ...
%         bar(bin_centers, neuron_params(n).h_norm);
end

save init_params pd gamma_const beta_const db training_period

% figure, plot(angle(:,1),'r')
list = [(1:158)' gamma_const' beta_const'];
disp(list);