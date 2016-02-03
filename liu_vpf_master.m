function outputData = liu_vpf_master(inputData,liuParams)

global nNeurons blocksz nParticles pNames hh

blocksz    = liuParams.blocksz;
nParticles = liuParams.nParticles;

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
end
rpfx.posn(:,1) = repmat(posn(1),nParticles,1);

% assign initial weights
g = (1/nParticles)*ones(nParticles,1);

% initialize theta particle filter
xEstimate.mu    = zeros(gen,nNeurons);
xEstimate.posn  = zeros(gen,1);
for i = 1:nNeurons
    xEstimate.mu(1,i)    = mean(rpfx.mu{i}(:,1));
end
xEstimate.posn(1) = mean(rpfx.posn(:,1));

pNames = {'mu','posn'};

for blockIndex = 1:gen-1
    fprintf('starting blockIndex %3i of %3i\n',blockIndex+1,gen);

    clear x
    for i = 1:nNeurons
        x.mu(:,i)  = rpfx.mu{i}(:,blockIndex);
        x.alpha(i) = 3.5; % changed from 4.0
        x.sigma(i) = 15;
    end
    x.posn = rpfx.posn(:,blockIndex);

    [sMu,prior] = liu_vpf_step1(x);
    g           = liu_vpf_step2a(x,N,g,blockIndex,sMu,dt,prior);

    % use the weights to resample
    [sMu,x,rpfx] = liu_vpf_step2b(g,sMu,x,rpfx);
     xNew.mu     = sMu.mu;
     xNew.posn   = sMu.posn;
    [rpfx]       = liu_assign_next_particle(rpfx,xNew,blockIndex);
    [xEstimate]  = liu_assign_next_estimate(xEstimate,g,rpfx,blockIndex);

end

outputData.rpfx = rpfx;
outputData.xEstimate = xEstimate;

