function outputData = liu_master(inputData,liuParams)

global nNeurons blocksz nParticles pNames hh

blocksz    = liuParams.blocksz;
nParticles = liuParams.nParticles;
nRepeats   = liuParams.nRepeats;

x_true   = inputData.x;
dt       = inputData.dt;
posn     = inputData.posn;
N        = inputData.N;
nNeurons = inputData.nNeurons;

gen  = floor(length(posn)/blocksz);
hh   = .01; % a constant of our choosing

% memory allocation
for i = 1:nNeurons
    rpfx.mu{i}    = zeros(nParticles,gen);
end
rpfx.posn = zeros(nParticles,gen);

% initialization
for i = 1:nNeurons
        rpfx.mu{i}(:,1)    = x_true.mu(i,1)*ones(nParticles,1); % our initial guess for mu
        %     rpfx.mu{i}(:,1)    = 150*ones(nParticles,1); % our initial guess for mu
        %     rpfx.mu{i}(:,1)    = 150*rand(nParticles,1); % our initial guess for mu
end
rpfx.posn(:,1) = repmat(posn(1),nParticles,1);
% rpfx.posn(:,1) = repmat(150,nParticles,1);
% rpfx.posn(:,1) = 150*rand(nParticles,1);

% assign initial weights
w = (1/nParticles)*ones(nParticles,1);

% initialize theta particle filter
xEstimate.mu    = zeros(gen,nNeurons);
xEstimate.posn  = zeros(gen,1);
for i = 1:nNeurons
    xEstimate.mu(1,i)    = mean(rpfx.mu{i}(:,1));
end
xEstimate.posn(1) = mean(rpfx.posn(:,1));

pNames = {'mu','posn'};
% h = [];
% close all; figure('windowstyle','docked');

for blockIndex = 1:gen-1
    fprintf('starting blockIndex %3i of %3i\n',blockIndex+1,gen);

%     subBlock_start = (blockIndex-1)*blocksz+1;
%     subBlock_end = (blockIndex)*blocksz;

    clear x
    for i = 1:nNeurons
        x.mu(:,i)  = rpfx.mu{i}(:,blockIndex);
%         x.alpha(i) = 4;
        x.alpha(i) = 3.5;
        x.sigma(i) = 15;
        %                 x.alpha(i) = mean(x_true.alpha(i,subBlock_start:subBlock_end));
        %                 x.sigma(i) = mean(x_true.sigma(i,subBlock_start:subBlock_end));
    end
    x.posn = rpfx.posn(:,blockIndex);

    [sMu,sTheta,m]=liu_step1(x);

    for qq = 1:nRepeats % this iterates steps 2-5 (liu) so you get better randomness

        g = liu_step2a(x,N,w,blockIndex,sMu,dt);

        % use the weights to resample
        [sTheta,sMu,x,rpfx,m] = liu_step2b(g,sTheta,sMu,x,rpfx,m);

        [sThetaNew] = liu_step3(sTheta,m);
        [xNew]      = liu_step4(sThetaNew,x);

        [rpfx] = liu_assign_next_particle(rpfx,xNew,blockIndex);

        w = liu_step5(N,blockIndex,x,xNew,sMu,dt);

    end

    [xEstimate] = liu_assign_next_estimate(xEstimate,w,rpfx,blockIndex);
    
%     delete(h);
%     currTime = blockIndex*blocksz*dt;
%     tt = (1:length(xEstimate.posn)) * dt * blocksz;
%     h = plot(t(t<currTime),posn(t<currTime),tt,xEstimate.posn,'-o');
%     xlim([0 tmax]); ylim([0 300]);
%     pause(0.1);
end

outputData.rpfx = rpfx;
outputData.xEstimate = xEstimate;
% %%
% close all; figure('windowstyle','docked');
% t = (1:length(posn))*dt;
% tt = (1:length(outputData.xEstimate.posn)) * dt * blocksz;
% plot(t,posn,tt,outputData.xEstimate.posn,'-o');
% % hold on; plot(tt,outputData.xEstimate.posn,'-ro');
% xlim([0 maxGen*blocksz*dt]);

