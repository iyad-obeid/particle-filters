clc
clear all
close all
load input_data

nParticles = 100;
blockSize  = 20;
nIterations= floor(length(t)/blockSize);

var_mu1 = 0.0001;
var_mu2 = 0.00001;

rho_const = rho_start;
sigma_squared_const = sigma_squared_start;
s_ds = downsample(s, blockSize);

P(1:nParticles) = struct('s',s_ds,...
    'rho', rho_const*ones(nIterations, nNeurons),...
    'mu',zeros(nIterations, nNeurons),...
    'sigma_squared', sigma_squared_const*ones(nIterations, nNeurons),...
    'g',zeros(nIterations,1),...
    'likelihood',zeros(nIterations,1),...
    'resampled_likelihood',zeros(nIterations,1),...
    'w',ones(nIterations,1));
x_estimate.mu= zeros(nIterations, nNeurons);

% initialize first iteration
for p =1:nParticles
    P(p).mu(1,:)= mu(1,:);
    x_estimate.mu(1, :) = mu_start;
end

tic;
for k = 2:nIterations
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % Step1: draw random samples  % % % % % % %
    for p = 1:nParticles
        mu_estimate = P(p).mu(k-1,:)/pi;
        P(p).mu(k,:) = pi*random_beta(mu_estimate, var_mu1);
        temp.mu(p,:) = P(p).mu(k,:);
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
        P(p).mu(k,:)=temp.mu(resampleIndicies(p),:);
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % Ste 4: draw random samples  % % % % % % %
    for p = 1:nParticles
        mu_estimate = P(p).mu(k,:)/pi;
        P(p).mu(k,:) = pi*random_beta(mu_estimate, var_mu2);
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
        x_estimate.mu(k,:) = P(p).mu(k,:)*P(p).w(k,1) + x_estimate.mu(k,:);
    end
end

tElapsed = toc;
disp(['Elapsed time = ', num2str(tElapsed/60), 'minutes']);

t2 = 0:dt*blockSize:tmax-dt;
figure, plot(t2, x_estimate.mu(:,1), t , mu(:,1), t, s, 'LineWidth',2),...
    ylim([0 pi]), grid on;

save results