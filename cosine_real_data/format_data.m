clear;close all;clc;
load paco071508a.mat

fs_kinematics = 36e3;
fs_neurons = 1e3;
t = (1:length(AD33))/fs_neurons;
% kinematics(1,:) = AD33; % shoulder angle
% kinematics(2,:) = AD34; % elbow angle
% kinematics(3,:) = AD35;
% kinematics(4,:) = AD36;
% kinematics(5,:) = AD37;
% kinematics(6,:) = AD38;
% kinematics(7,:) = AD39;
% kinematics(8,:) = AD40;
% x = 0.1:0.1:(pi-0.1);

names = who;
p = 0;
for k = 1: length(names)
    if strncmp(names(k), 'sig', 3) == 1
        p = p + 1;
        s = char(strcat('temp = round(', names(k),'*fs_kinematics);'));
        eval(s);
        firings = false(1,((fs_kinematics/fs_neurons)*length(t)));
        firings(temp) = 1;
        firings_decimated = false(1, length(AD33));
        f = floor(find(firings)/36);
        firings_decimated(f) = 1;
        s = (['N', num2str(p), ' = firings;']);
        eval(s);
        eval(['N', num2str(p), '_resampled = firings_decimated;']);
%        eval(['h = find(N', num2str(p), '_resampled.*kinematics(1,:));']);
%         figure(p), set(gcf, 'WindowStyle', 'docked'), ...
%             hist(kinematics(1,h), x);
        save(['neuron', num2str(p)],...
            ['N', num2str(p)], ['N', num2str(p),'_resampled']);
        eval(['N_total(1,', num2str(p) ,') = sum(N', num2str(p),');']);
        s = ['clear N', num2str(p), ' N',num2str(p),'_resampled',';'];
        eval(s);
    end
end

kinematics(1,:) = AD33; % shoulder angle
kinematics(2,:) = AD34; % elbow angle
% kinematics(3,:) = AD35;
% kinematics(4,:) = AD36;
% kinematics(5,:) = AD37;
% kinematics(6,:) = AD38;
% kinematics(7,:) = AD39;
% kinematics(8,:) = AD40;

dt = 1/fs_neurons;
save('kinematics', 'kinematics', 'dt', 'N_total');

figure, plot(t, kinematics(1,:));
figure, plot(t, kinematics(2,:));
% figure, plot(t, kinematics(3,:));
% figure, plot(t, kinematics(4,:));
% figure, plot(t, kinematics(5,:));
% figure, plot(t, kinematics(6,:));
% figure, plot(t, kinematics(7,:));
% figure, plot(t, kinematics(8,:));
