clc
clear all
close all
load input_data

nParticles = 50;
blockSize  = 20;
nIterations= floor(length(t)/blockSize);

kappa.theta_p1 = 300;
kappa.theta_p2 = 300;

gamma_const = gamma_start;
beta_const = beta_start;
s_ds = downsample(s, blockSize);

% P(1:nParticles) = struct('s',s_ds,...
%     'gamma', gamma_const*ones(nIterations, nNeurons),...
%     'theta_p',zeros(nIterations, nNeurons),...
%     'beta', beta_const*ones(nIterations, nNeurons),...
%     'g',zeros(nIterations,1),...
%     'likelihood',zeros(nIterations,1),...
%     'resampled_likelihood',zeros(nIterations,1),...
%     'w',ones(nIterations,1));
% x_estimate.theta_p= zeros(nIterations, nNeurons);


P_new(1:nParticles) = struct('s',0,...
    'gamma', gamma_const*ones(1, nNeurons),...
    'theta_p',zeros(1, nNeurons),...
    'beta', beta_const*ones(1, nNeurons),...
    'g',zeros(1,1),...
    'likelihood',zeros(1,1),...
    'resampled_likelihood',zeros(1,1),...
    'w',ones(1,1));
P_old = P_new;

x_estimate.theta_p= zeros(nIterations, nNeurons);


% initialize first iteration
for p =1:nParticles
    P(p).theta_p(1,:)= theta_p(1,:);
    x_estimate.theta_p(1, :) = theta_p_start;
end

tic;
for k = 2:nIterations
    disp(['block number ', num2str(k), ' of ', num2str(nIterations)]);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % Step1: draw random samples  % % % % % % %
    for p = 1:nParticles
        for n = 1:nNeurons
          P(p).theta_p(k,n)   = mod(randraw('vonmises', [P(p).theta_p(k-1,n), kappa.theta_p1], 1), 2*pi);
        end 
        temp.theta_p(p,:) = P(p).theta_p(k,:);
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % Step2: compute likelihoods % % % % % % % %
    norm_sum = 0;
    for p = 1:nParticles
        particle = P(p);
        blockStart = (k-1)*blockSize + 1;
        blockEnd   = blockStart + blockSize - 1;
        nFirings = sum(N(blockStart:blockEnd,:),1);
        P(p).likelihood(k,1)=compute_likelihood(particle, nFirings, blockSize, k, dt);
        P(p).g(k,1) = P(p).w((k-1),1)*P(p).likelihood(k,1);
        norm_sum = norm_sum + P(p).g(k,1);
    end
    % normalize
    prob_vector = zeros(1, nParticles);
    for p = 1:nParticles
        P(p).g(k,1) = P(p).g(k,1)/norm_sum;
        prob_vector(1,p) = P(p).g(k,1);
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % Step 3: resample % % % % % % % % % % % % %
    resampleIndicies = randsample((1:nParticles), nParticles, 'true', prob_vector);
    for p = 1:nParticles
        P(p).theta_p(k,:)=temp.theta_p(resampleIndicies(p),:);
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % Ste 4: draw random samples  % % % % % % %
    for p = 1:nParticles
        for n = 1:nNeurons
            P(p).theta_p(k,n) = mod(randraw('vonmises', [P(p).theta_p(k,n), kappa.theta_p2], 1), 2*pi);
        end
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % Step 5: Compute resampled likelihoods % %
    norm_sum = 0;
    for p = 1:nParticles
        particle = P(p);
        P(p).resampled_likelihood(k,1) = compute_likelihood(particle, nFirings, blockSize, k, dt);
        P(p).w(k,1) = P(p).resampled_likelihood(k,1)/P(p).likelihood(k,1);
        norm_sum = norm_sum + P(p).w(k,1);
    end
    % normalize
    for p = 1:nParticles
        P(p).w(k,1) = P(p).w(k,1)/norm_sum;
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % Step 6: Compute weighted estimate % % % % % % % % % % %
    for p = 1:nParticles
        x_estimate.theta_p(k,:) = P(p).theta_p(k,:)*P(p).w(k,1) + x_estimate.theta_p(k,:);
    end
end

tElapsed = toc;
disp(['Elapsed time = ', num2str(tElapsed/60), 'minutes']);
t2 = 0:dt*blockSize:tmax-dt;
figure, plot(t2, x_estimate.theta_p(:,1), t , theta_p(:,1), t, s, 'LineWidth',2),...
 %   ylim([0 2*pi]), grid on;
save results