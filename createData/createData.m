clear; clf;

dt = 20e-3;
t = dt:dt:500;
posn = 150*sawtooth(2*pi/(5)*t,.5)+150;

alpha=log(30);
sigma = 11;

x1 = 150; % in real life, this is mu, the optimal firing posn of the neuron

lambda=exp(alpha-((posn-x1).^2)./(2*sigma.^2)); %conditional intesity function

N1 = poissrnd(lambda*dt);
N1 = N1>0;

subplot(2,1,1);
plot(t,posn);
hold on
r = find(N1);
plot(t(r),posn(r),'ro');

%%

x2 = zeros(size(t));
r1 = t<250;
x2(r1) = polyval([(450-100)/250 100],t(r1));
r2 = t>=250;
x2(r2) = polyval([(100-450)/250 800],t(r2));

lambda = exp(alpha - (posn-x2).^2 / (2*sigma^2));

N2 = poissrnd(lambda*dt);
N2 = N2>0;

subplot(2,1,2);
plot(t,posn);
hold on
r = find(N2);
plot(t(r),posn(r),'ro');


%%

clear; clf;

dt = 20e-3;
t = dt:dt:500;
posn = 150*sawtooth(2*pi/(5)*t,.5)+150;

x.alpha=linspace(log(30),log(10),length(t));
x.sigma=linspace(9,13,length(t));
x.mu = linspace(50,250,length(t));

lambda = exp(x.alpha - (posn-x.mu).^2 ./ (2*x.sigma.^2));

N = poissrnd(lambda*dt);
N = N>0;

plot(t,posn);
hold on
r = find(N);
plot(t(r),posn(r),'ro');

% clear r;
% save data_feb27

%%
clear;clf;
dt = 20e-3;
t = dt:dt:500;
posn = 150*sawtooth(2*pi/(5)*t,.5)+150;

x.alpha = log(30) * ones(size(posn));
x.sigma = 11 * ones(size(posn));

x.mu = zeros(size(t));
r1 = t<250;
x.mu(r1) = polyval([(450-100)/250 100],t(r1));
r2 = t>=250;
x.mu(r2) = polyval([(100-450)/250 800],t(r2));

lambda = exp(x.alpha - (posn-x.mu).^2 ./ (2*x.sigma.^2));

N = poissrnd(lambda*dt);
N = N>0;

plot(t,posn);
hold on
r = find(N);
plot(t(r),posn(r),'ro');

clear r;
save data_mar04_v1

%%
clear;clf;
dt = 2e-3;
t = dt:dt:30;
posn = 150*sin(2*pi*t/5) + 150;
nNeurons = 50;
x.alpha = log(rand(nNeurons,1)*50) * ones(size(posn));
x.sigma = (rand(nNeurons,1)*10 + 10) * ones(size(posn));
x.mu = zeros(size(x.alpha));
for i = 1:nNeurons
   kg = 1;
   while (kg)
       mStart = rand * 300;
       mEnd = rand * 300;
       if abs(mStart-mEnd)<50
           kg = 0;
       end
   end
   x.mu(i,:) = linspace(mStart,mEnd,length(posn));
end

ptmp = ones(nNeurons,1) * posn;
lambda = exp(x.alpha - (ptmp-x.mu).^2 ./ (2*x.sigma.^2));

N = poissrnd(lambda*dt);
N = N>0;

inputData.x = x;
inputData.dt = dt;
inputData.posn = posn;
inputData.N = N;
inputData.nNeurons = nNeurons;

liuParams.blocksz = 25;
liuParams.nParticles = 100;

save ../data/data_jan06_v2 inputData liuParams
