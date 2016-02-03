function inputData = create_distractor_model(nNeurons,pctPrimary,tMax)

dt         = 2e-3;
t          = dt:dt:tMax;
model      = @(a,m,s,x) exp(a- ((x-m).^2)./(2*s.^2) );

max_firing = 100; %spikes per second

posn       = 150*sin(2*pi/5*t) + 150;

% alpha and sigma for each neuron are constant in time
% mu changes linearly in time with a max slope of +/- 50 units per 100
% seconds

alpha = log(max_firing)*rand(nNeurons,1); %unif distro bt 0 & max_firing
alpha = repmat(alpha,1,length(t));

sigma = 10*rand(nNeurons,1)+10; %uniformly distro bt 10 & 20
sigma = repmat(sigma,1,length(t));

mu_slope   = (2*rand(nNeurons,1)-1)/2; %unif distro bt -0.5 & 0.5
mu_initial = 400*rand(nNeurons,1)-50; %unif distro bt -50 & 350
mu_final   = (tMax*mu_slope)+mu_initial;
mu         = zeros(nNeurons,length(t));

for L = 1:nNeurons
    mu(L,:) = linspace(mu_initial(L,1),mu_final(L,1),length(t));
end

nPrimary = floor(nNeurons * pctPrimary);
r1 = 1:nPrimary;
r2 = nPrimary + 1 : nNeurons;
lambda = zeros(nNeurons,length(t));

lambda(r1,:)    = model(alpha(r1,:),mu(r1,:),sigma(r1,:),repmat(posn,nPrimary,1) );
for i = r2
    lambda(i,:) = model(alpha(i,:),mu(i,:),sigma(i,:),create_distractor(length(t)));
end
N = poissrnd(lambda*dt);
N = N>0;

inputData.dt              = dt;
inputData.t               = t;
inputData.posn            = posn;
inputData.N               = N;
inputData.x.alpha         = alpha;
inputData.x.mu            = mu;
inputData.x.sigma         = sigma;
inputData.nNeurons        = nNeurons;

function distractor = create_distractor(nPts)
dx         = 1.5;
flag        =  0;
while (flag == 0)
    distractor = 150 + cumsum(dx*randn(1,nPts));
    if (max(distractor)<=300 && min(distractor)>=0)
        flag = 1;
    end
end


