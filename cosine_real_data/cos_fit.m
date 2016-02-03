clear; clf; close all;
load paco071508a;
t = (1:length(AD33))/1e3;
bin_edges = linspace(-2*pi,0,50)';

AD33(AD33>-0.5) = mean(AD33);

s = whos('sig*');
pdVec = zeros(1,length(s));

for i = 1:length(s)
    figure('windowstyle','docked');
    disp(i);
    sp_times = eval(s(i).name);
    ad33 = interp1(t,AD33,sp_times);
    cnt1 = histc(ad33,bin_edges);
    cnt2 = histc(AD33,bin_edges);
    cnt = cnt1./cnt2;
    bar(bin_edges,cnt,'histc');
    
    bin = bin_edges(~isnan(cnt));
    cnt = cnt(~isnan(cnt));
    
    [~,j] = max(cnt);
    pd = bin(j);
    A = [cos(bin-pd) ones(size(bin))]\cnt;
    y = A(1)*cos(bin-pd) + A(2);
    hold on
    plot(bin,y,'r','linewidth',2);
    pdVec(i) = pd;
%     drawnow;
%     pause(0.1);
end

