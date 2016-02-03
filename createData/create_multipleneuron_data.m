function inputData = create_multipleneuron_data(nNeurons,T,varargin)

FLAT_MU = 0;
if nargin>2
    for i = 3:nargin
        ind = i -2;
        if strcmp(varargin{ind},'flat')
            FLAT_MU = 1;
        end
    end
end

max_firing = 100; %spikes per second
dt = 2e-3; %seconds

t = dt:dt:T; %seconds

% alpha and sigma for each neuron are constant in time
% mu changes linearly in time with a max slope of +/- 50 units per 100
% seconds
alpha=log(max_firing)*rand(nNeurons,1); %unif distro bt 0 & max_firing
alpha=repmat(alpha,1,length(t));
sigma=10*rand(nNeurons,1)+10; %uniformly distro bt 10 & 20
sigma=repmat(sigma,1,length(t));

mu_slope=(2*rand(nNeurons,1)-1)/2; %unif distro bt -0.5 & 0.5
mu_initial=400*rand(nNeurons,1)-50; %unif distro bt -50 & 350
mu_final=(T*mu_slope)+mu_initial;
mu=zeros(nNeurons,length(t));

if (FLAT_MU==0)
    for L=1:nNeurons
        mu(L,:)=linspace(mu_initial(L,1),mu_final(L,1),length(t));
    end
else
    mu=repmat(mu_initial,1,length(t));
end
model = @(a,m,s,x) exp(a- ((x-m).^2)./(2*s.^2) );

posn = 150*sin(2*pi/5 * t) + 150;
lambda = model(alpha,mu,sigma,repmat(posn,nNeurons,1));
N = poissrnd(lambda*dt);
N = N>0;

inputData.dt = dt;
inputData.t = t;
inputData.posn = posn;
inputData.N = N;
inputData.x.alpha = alpha;
inputData.x.mu = mu;
inputData.x.sigma = sigma;
inputData.nNeurons = nNeurons;

% save data02 inputData;
