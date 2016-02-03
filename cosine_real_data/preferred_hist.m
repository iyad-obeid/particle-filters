clc
clear
load kinematics

kinematics = kinematics(:, 2);
kinematics((kinematics>pi)) = mean(kinematics);
N = false(length(kinematics), 158);

for n = 1:158
    eval(['load neuron', num2str(n)]);
    eval(['temp = N', num2str(n),'(1:length(kinematics));']);
    N(:, n) = temp';
    %eval(['N(:,', num2str(n),') = N', num2str(n), '(1:length(kinematics(:,2)));']);
    eval(['clear N', num2str(n), ' N', num2str(n)]);
end
F_elbow = zeros(size(N));
for n = 1:158
    F_elbow(:,n) = N(:,n).*kinematics;
end
save F_elbow F_elbow
clear N N_total
db = 0.1;
bin_centers = db:db:pi;

kin_hist = hist(kinematics, bin_centers);

neuron_params(1:158)= struct('pos',0);
pd = zeros(1, 158);
gamma_const = zeros(1, 158);
beta_const = zeros(1, 158);

for n = 1:158
    neuron_params(n).pos = find(F_elbow(:,n));
    neuron_params(n).h = hist(kinematics(neuron_params(n).pos), bin_centers);
    neuron_params(n).h_norm = (neuron_params(n).h)./kin_hist;
    
    [~, index] = max(neuron_params(n).h_norm);
    pd(n) = bin_centers(index);
    neuron_params(n).max_rate = (neuron_params(n).h(index))/(kin_hist(index)*dt);
    [~, index] = min(neuron_params(n).h_norm);
    neuron_params(n).min_rate = (neuron_params(n).h(index))/((kin_hist(index)*dt));
    
    gamma_const(1, n) = (neuron_params(n).max_rate - neuron_params(n).min_rate)/2;
    beta_const(1, n) = gamma_const(1, n) + neuron_params(n).min_rate;
%     
%     figure(n), set(gcf,'WindowStyle','docked'), ...
%         bar(bin_centers, neuron_params(n).h_norm);
end

save preferred_hist pd gamma_const beta_const db
%clear F_elbow
% figure, plot(kinematics(:,1),'r')
list = [(1:158)' gamma_const' beta_const'];
disp(list);

figure, hist(pd, bin_centers)