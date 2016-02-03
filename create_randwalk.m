% This file creates the random walk and then runs the BAPF

clear;
dt = 2e-3;
tmax = 30;

dx = 1.5;
nNeurons = 10;

t = dt:dt:tmax;
nPts = length(t);
startPosn = 150;
minPosn   = 0;
maxPosn   = 300;

% Define mouse trajectory
flag = 1;
while (flag)
    posn =  startPosn + cumsum(dx*randn(size(t)));
    if ((max(posn)<maxPosn) && (min(posn)>minPosn))
        flag = 0;
    end
end

% Define Alphas
x.alpha = zeros(nNeurons,nPts);
for i = 1:nNeurons
    flag = 1;
    while (flag)
        aStart = 50*rand;
        aEnd = 50*rand;
        if abs(aStart-aEnd)<5
            flag = 0;
        end
    end
    x.alpha(i,:) = log(linspace(aStart,aEnd,nPts));
end

% Define Sigmas
x.sigma = zeros(nNeurons,nPts);
for i = 1:nNeurons
    flag = 1;
    while (flag)
        sStart = 10 + 10*rand;
        sEnd   = 10 + 10*rand;
        if abs(sStart-sEnd)<5
            flag = 0;
        end
    end
    x.sigma(i,:) = linspace(sStart,sEnd,nPts);
end

% Define Mus
x.mu = zeros(nNeurons,nPts);
for i = 1:nNeurons
    mStart = -50 + 400*rand;
    mSlope = (2*rand-1)/2;
    x.mu(i,:) = mStart + mSlope*t;
end

% Define N
posnMatrix = ones(nNeurons,1) * posn;
lambda = exp(x.alpha - ((x.mu - posnMatrix).^2) ./ (2*x.sigma.^2));
N = poissrnd(lambda*dt);

% Define constants and package them for easy transport to "liu_master"
liuParams.blocksz = 25;
liuParams.nParticles = 100;
liuParams.nRepeats = 5;

inputData.x = x;
inputData.dt       = dt;
inputData.posn     = posn;
inputData.N = N;
inputData.nNeurons = nNeurons;

% run the BAPF simulation
outputData = liu_master(inputData,liuParams);

% capture output and create plot
posnEstimate = outputData.xEstimate.posn;
dt2 = liuParams.blocksz * dt;
tEstimate = dt2 : dt2 : tmax;
clf;
plot(t,posn,tEstimate,posnEstimate); axis([0 30 0 300]);
