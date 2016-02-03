clc
clear all
load kinematics
load initial_estimates

% There are 687974 timesteps [t=0-687.9 seconds] in the original data. We can elect to work
% with only a subset if we choose

decode_period   = (300:dt:(330-dt));
tStart          =  decode_period(1);
tEnd            = decode_period(end);
t               = (tStart:dt:tEnd-dt)'; % almost the same as decode_period
startSample     = (decode_period(1))/dt;
endSample       = (decode_period(end))/dt;
samples         = (startSample:endSample)';
s               = kinematics(samples, 2);

clear kinematics blockSize



% There are 158 neurons altogether. 1-79 come from one electrode bundle and
% 80-158 come from another electrode bundle. It won't make sense to try to
% mix them together. Pick one or the other!
% neuron_list = 1:79;
neuron_list     = 80:158;
nNeurons        = length(neuron_list);
load neuronData % This is where the spike data are saved
N               = N(samples,neuron_list); % just keep the spikes we'll be needing

blockSize   = 100;
s_ds        = downsample(s, blockSize);
nIterations = floor(length(s)/blockSize);
nParticles  = 1000;
sd.s1       = 0.15;
sd.s2       = 0.05;
sd.theta_p1 = 0.001;
sd.theta_p2 = 0.001;
sd.gamma1   = 0.01;
sd.gamma2   = 0.001;
sd.beta1    = 0.01;
sd.beta2    = 0.001;

P_new(1:nParticles) = struct(...
    's'                     ,0,...
    'gamma'                 ,zeros(1, nNeurons),...
    'theta_p'               ,zeros(1, nNeurons),...
    'beta'                  ,zeros(1, nNeurons),...
    'g'                     ,0,...
    'likelihood'            ,0,...
    'resampled_likelihood'  ,0,...
    'w'                     ,1);
P_old = P_new;

% initialize first iteration
s_init          = s(1);
theta_p_init    = theta_p_final(neuron_list);
gamma_init      = gamma_final(neuron_list);
beta_init       = beta_final(neuron_list);
for p =1:nParticles
    P_old(p).theta_p(1,:)   = theta_p_init;
    P_old(p).gamma(1,:)     = gamma_init;
    P_old(p).beta(1,:)      = beta_init;
    P_old(p).s              = s_init;
    P_new(p).gamma  = P_old(p).gamma;
    P_new(p).beta   = P_old(p).beta;
end
x_estimate.theta_p  = zeros(nIterations, nNeurons);
x_estimate.gamma    = zeros(nIterations, nNeurons);
x_estimate.beta     = zeros(nIterations, nNeurons);
x_estimate.s        = zeros(nIterations, 1       );

tic;
for k = 1:nIterations
    disp(['block number ', num2str(k), ' of ', num2str(nIterations)]);
 
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % Step1: draw random samples  % % % % % % %
    for p = 1:nParticles
        
        flag = true;
        while flag
            v = P_old(p).s + sd.s1*randn;
            if v > 0 && v < pi
                P_new(p).s = v;
                flag = false;
            end
        end
        
        flag = true;
        while flag
            v = P_old(p).theta_p + sd.theta_p1*randn(1,nNeurons);
            if all(v > 0 & v < pi)
                P_new(p).theta_p = v;
                flag = false;
            end
        end
        
        temp.s(p)           = P_new(p).s;
        temp.theta_p(p,:)   = P_new(p).theta_p;

    end
    
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % Step2: compute likelihoods % % % % % % % %
    norm_sum = 0;
    for p = 1:nParticles
        particle    = P_new(p);
        blockStart  = (k-1)*blockSize + 1;
        blockEnd    = blockStart + blockSize - 1;
        nFirings    = sum(N(blockStart:blockEnd,:),1);
        
        P_new(p).likelihood = compute_likelihood(particle, nFirings, blockSize, dt);
        P_new(p).g          = P_old(p).w * P_new(p).likelihood;
        norm_sum = norm_sum + P_new(p).g;
    end
    
    % normalize
    prob_vector = zeros(1, nParticles);
    for p = 1:nParticles
        P_new(p).g      = P_new(p).g/norm_sum;
        prob_vector(p)  = P_new(p).g;
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % Step 3: resample % % % % % % % % % % % % %
    resampleIndicies = randsample(nParticles, nParticles, true, prob_vector);
    for p = 1:nParticles
        P_new(p).theta_p    = temp.theta_p(resampleIndicies(p),:);
        P_new(p).s          = temp.s      (resampleIndicies(p)  );
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % Step 4: draw random samples  % % % % % % %
    for p = 1:nParticles
        flag = true;
        while flag
            v = P_new(p).s(1,1)+ sd.s2*randn;
            if v > 0 && v < pi
                P_new(p).s(1,1) = v;
                flag = false;
            end
        end
        
        flag = true;
        while flag
            v = P_new(p).theta_p(1,:)+ sd.theta_p2*randn(1,nNeurons);
            if all(v > 0 & v < pi)
                P_new(p).theta_p(1,:) = v;
                flag = false;
            end
        end
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % Step 5: Compute resampled likelihoods % %
    norm_sum = 0;
    for p = 1:nParticles
        particle = P_new(p);
        P_new(p).resampled_likelihood   = compute_likelihood(particle, nFirings, blockSize, dt);
        P_new(p).w                      = P_new(p).resampled_likelihood / P_new(p).likelihood;
        norm_sum                        = norm_sum + P_new(p).w;
    end
    % normalize
    for p = 1:nParticles
        P_new(p).w = P_new(p).w/norm_sum;
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % Step 6: Compute weighted estimate % % % % % % % % % % %
    for p = 1:nParticles
        x_estimate.s(k)         = P_new(p).s       * P_new(p).w + x_estimate.s      (k  );
        x_estimate.theta_p(k,:) = P_new(p).theta_p * P_new(p).w + x_estimate.theta_p(k,:);
        x_estimate.gamma  (k,:) = P_new(p).gamma   * P_new(p).w + x_estimate.gamma  (k,:);
        x_estimate.beta(k,:)    = P_new(p).beta    * P_new(p).w + x_estimate.beta   (k,:);
    end
    P_old = P_new;
end
toc
%tElapsed = toc;
%disp(['Elapsed time = ', num2str(tElapsed/60), 'minutes']);

t2 = tStart:dt*blockSize:tEnd-dt;

p = polyfit(s_ds, x_estimate.s(:,1),1);
f = polyval(p, s_ds);
sserr = sum((x_estimate.s(:,1) - f).^2);
sstot = sum((x_estimate.s(:,1) - mean(x_estimate.s(:,1))).^2);
r2 = 1-sserr/sstot;
disp(sqrt(r2));
%x_estimate.s(201:220) = x_estimate.s(201:220)*.85;
L1 = ['# neurons = ', num2str(nNeurons), ...
    ' # particles = ', num2str(nParticles), ...
    ' block size = ', num2str(blockSize)];
L2 = ['\sigma_{s1} = ' ,num2str(sd.s1), '   \sigma_{s2} = ', num2str(sd.s2)];
L3 = ['r = ', num2str(sqrt(r2))];

figure, plot(t2, x_estimate.s,'k-.','LineWidth',1), grid on, hold on, ...
    plot(t2, s_ds, 'k','LineWidth',2), ...
    xlabel('time (s)'), ylabel('\phi (rad)'), ...
    xlim([tStart tEnd]), ylim([0 pi]), ...
    legend('\phi_{est}', '\phi_{true}','Location','SouthEast'),
    %title({L1; L2; L3});

mse_BAPF = mean((x_estimate.s-s_ds).^2);
% save results