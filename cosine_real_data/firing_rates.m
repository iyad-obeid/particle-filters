clear;close all;clc;
load paco071508a.mat

names = who;
p = 0;
for k = 1: length(names)
    if strncmp(names(k), 'sig', 3) == 1
        p = p + 1;
        s = char(strcat('temp = (', names(k),'*1000);'));
        eval(s);
        temp = temp/1000;
        d = diff(temp);
        r = 1./d(:,1);
        mean_fr(p) = mean(r);
       % figure(p), set(gcf, 'WindowStyle','docked'), plot(r);
    end
end