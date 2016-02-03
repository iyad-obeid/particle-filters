function create_mixture_model(fname,nNeurons,pctPrimary,tMax)
% NOT COMPLETE - NEEDS EDITING!

dt         = 2e-3;
t          = 0:dt:tMax;

max_firing = 100; %spikes per second

driving1   = 150*sin(2*pi/5*t) + 150;
dx         = 1.5;
flag       = 0;
while (flag == 0)
    driving2 = 150 + cumsum(dx*randn(1,length(t)));
    if (max(driving2)<=300 && min(driving2)>=0)
        flag = 1;
    end
end

figure('windowstyle','docked');
plot(t,driving1,t,driving2)

w1    = 0.80;
w2    = 1 - w1;

model = @(a,m,s,x) exp(a- ((x-m).^2)./(2*s.^2) );

% alpha and sigma for each neuron are constant in time
% mu changes linearly in time with a max slope of +/- 50 units per 100
% seconds

alpha = log(max_firing)*rand(nNeurons*2,1); %unif distro bt 0 & max_firing
alpha = repmat(alpha,1,length(t));

sigma = 10*rand(nNeurons*2,1)+10; %uniformly distro bt 10 & 20
sigma = repmat(sigma,1,length(t));

mu_slope   = (2*rand(nNeurons*2,1)-1)/2; %unif distro bt -0.5 & 0.5
mu_initial = 400*rand(nNeurons*2,1)-50; %unif distro bt -50 & 350
mu_final   = (tMax*mu_slope)+mu_initial;
mu         = zeros(nNeurons*2,length(t));

for L = 1:2*nNeurons
    mu(L,:) = linspace(mu_initial(L,1),mu_final(L,1),length(t));
end

lam1   = model(alpha(1:nNeurons,:) ,mu(1:nNeurons,:), sigma(1:nNeurons,:), repmat(driving1,nNeurons,1));
lam2   = model(alpha(nNeurons+1:2*nNeurons,:) ,mu(nNeurons+1:2*nNeurons,:), sigma(nNeurons+1:2*nNeurons,:), repmat(driving2,nNeurons,1));
lambda = w1*lam1 + w2*lam2;

N = poissrnd(lambda*dt);
N = N>0;

inputData.dt              = dt;
inputData.t               = t;
inputData.posn            = driving1;
inputData.posn_distractor = driving2;
inputData.weights         = [w1 w2];
inputData.N               = N;
inputData.x.alpha         = alpha;
inputData.x.mu            = mu;
inputData.x.sigma         = sigma;
inputData.nNeurons        = nNeurons;

save(fname,inputData);


