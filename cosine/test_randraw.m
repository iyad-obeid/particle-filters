clc
clear

% y = randraw('vonmises', [pi, 5], [1 1e5]);
% 
% m = mean(y);
% v = var(y);
% 
% figure, hist(y, 200)



x = mod(randraw('vonmises', [2*pi, 200], 10000), 2*pi);