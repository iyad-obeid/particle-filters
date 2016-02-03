clear all ; clc;
% close all
global nNeurons blocksz nParticles pNames hh

blocksz    = 25;
nParticles = 100;
nRepeats   = 5;
load('data04');

x_true   = inputData.x;
dt       = inputData.dt;
posn     = inputData.posn;
N        = inputData.N;
nNeurons = inputData.nNeurons;

gen  = floor(length(posn)/blocksz);
hh   = .01; % a constant of our choosing
maxGen = 150;

% memory allocation
for i = 1:nNeurons
    rpfxA.mu{i}    = zeros(nParticles,gen);
    rpfxB.mu{i}    = zeros(nParticles,gen);
end
rpfxA.posn = zeros(nParticles,gen);
rpfxB.posn = zeros(nParticles,gen);

% initialization
for i = 1:nNeurons
    rpfxA.mu{i}(:,1)    = x_true.mu(i,1)*ones(nParticles,1); % our initial guess for mu
    rpfxB.mu{i}(:,1)    = x_true.mu(i,1)*ones(nParticles,1); % our initial guess for mu
    %     rpfx.mu{i}(:,1)    = 150*ones(nParticles,1); % our initial guess for mu
    %     rpfx.mu{i}(:,1)    = 150*rand(nParticles,1); % our initial guess for mu
end
rpfxA.posn(:,1) = repmat(posn(1),nParticles,1);
rpfxB.posn(:,1) = repmat(posn(1),nParticles,1);
% rpfx.posn(:,1) = repmat(150,nParticles,1);
% rpfx.posn(:,1) = 150*rand(nParticles,1);

% assign initial weights
wA = (1/nParticles)*ones(nParticles,1);
wB = (1/nParticles)*ones(nParticles,1);

% initialize theta particle filter
xEstimateA.mu    = zeros(gen,nNeurons);
xEstimateA.posn  = zeros(gen,1);
xEstimateB.mu    = zeros(gen,nNeurons);
xEstimateB.posn  = zeros(gen,1);
for i = 1:nNeurons
    xEstimateA.mu(1,i)    = mean(rpfxA.mu{i}(:,1));
    xEstimateB.mu(1,i)    = mean(rpfxB.mu{i}(:,1));
end
xEstimateA.posn(1) = mean(rpfxA.posn(:,1));
xEstimateB.posn(1) = mean(rpfxB.posn(:,1));

pNames = {'mu','posn'};
h = [];close all; figure('windowstyle','docked');
tmax = maxGen*dt*blocksz;    t = (1:length(posn))*dt;
for blockIndex = 1:maxGen-1
    for interleave = {'A','B'}
        fprintf('starting blockIndex %3i%c of %3i\n',blockIndex+1,interleave{1},maxGen);

        if strcmp(interleave,'A')
            subBlock_start = (blockIndex-1)*blocksz+1;
            subBlock_end = (blockIndex)*blocksz;
            rpfxLocal = rpfxA;
            wLocal = wA;
            xEstimateLocal = xEstimateA;
        else % case B
            subBlock_start = (blockIndex-1)*blocksz+1 + (blocksz/2);
            subBlock_end = (blockIndex)*blocksz + (blocksz/2);
            rpfxLocal = rpfxB;
            wLocal = wB;
            xEstimateLocal = xEstimateB;
        end
        % make sure we don't fall off the end
        
        clear x
        for i = 1:nNeurons
            x.mu(:,i)  = rpfxLocal.mu{i}(:,blockIndex);
            x.alpha(i) = 4;
            x.sigma(i) = 15;
        end
                
        x.posn = rpfxLocal.posn(:,blockIndex);
        [sMu,sTheta,m]=liu_step1(x);

        for qq = 1:nRepeats % this iterates steps 2-5 (liu) so you get better randomness

            g = liu_step2a(x,N,wLocal,blockIndex,sMu,dt);

            % use the weights to resample
            [sTheta,sMu,x,rpfx,m] = liu_step2b(g,sTheta,sMu,x,rpfxLocal,m);

            [sThetaNew] = liu_step3(sTheta,m);
            [xNew]      = liu_step4(sThetaNew,x);

            [rpfxLocal] = liu_assign_next_particle(rpfxLocal,xNew,blockIndex);

            wLocal = liu_step5(N,blockIndex,x,xNew,sMu,dt);

        end

        [xEstimateLocal] = liu_assign_next_estimate(xEstimateLocal,wLocal,rpfxLocal,blockIndex);

        if strcmp(interleave,'A')
           xEstimateA = xEstimateLocal;
           rpfxA = rpfxLocal;
           wA = wLocal;
        else
           xEstimateB = xEstimateLocal;
           rpfxB = rpfxLocal;
           wB = wLocal;
        end
    end
end


outputData.rpfxA = rpfxA;
outputData.rpfxB = rpfxB;
outputData.xEstimateA = xEstimateA;
outputData.xEstimateB = xEstimateB;


%%
close all; figure('windowstyle','docked');
t = (1:length(posn))*dt;
tt = (1:length(xEstimateA.posn))*dt*blocksz;
tt2 = (1:length(xEstimateA.posn))*dt*blocksz + (dt*blocksz/2);
plot(tt,xEstimateA.posn,'ro');
hold on; 
plot(tt2,xEstimateB.posn,'ko');
plot(t,posn,'g'); 
xlim([0 8]);
xlabel('time/s');
bigText;
