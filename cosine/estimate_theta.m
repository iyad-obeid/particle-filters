clc
clear all
close all
load input_data

nParticles = 200;
blockSize  = 30;
sd.theta_p1 = 0.01;
sd.theta_p2 = 0.01;
nIterations= floor(length(t)/blockSize);

gamma_const = 8;
beta_const = gamma_const;
% theta_p_init = theta_p_start;
theta_p_init = 2.5;

P_new(1:nParticles) = struct('s',0,...
    'gamma', gamma_const*ones(1, nNeurons),...
    'theta_p',zeros(1, nNeurons),...
    'beta', beta_const*ones(1, nNeurons),...
    'g',zeros(1,1),...
    'likelihood',zeros(1,1),...
    'resampled_likelihood',zeros(1,1),...
    'w',ones(1,1));
P_old = P_new;
s_ds = downsample(s, blockSize);
x_estimate.theta_p= zeros(nIterations, nNeurons);

% initialize first iteration
for p =1:nParticles
    P_old(p).theta_p(1,:)= theta_p_init;
    x_estimate.theta_p(1, :) = theta_p_init;
end

tic;
for k = 2:nIterations
    disp(['block number ', num2str(k), ' of ', num2str(nIterations)]);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % Step1: draw random samples  % % % % % % %
    
    for p = 1:nParticles
        P_new(p).s = s_ds(k);
        flag = true;
        while flag
            v = P_old(p).theta_p(1,:)+ sd.theta_p1*randn(1, nNeurons);
            if v > 0 & v < pi
                P_new(p).theta_p(1,:) = v;
                flag = false;
            end
        end
        
        temp.theta_p(p,:) = P_new(p).theta_p(1,:);
        
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % Step2: compute likelihoods % % % % % % % %
    norm_sum = 0;
    for p = 1:nParticles
        particle = P_new(p);
        blockStart = (k-1)*blockSize + 1;
        blockEnd   = blockStart + blockSize - 1;
        nFirings = sum(N(blockStart:blockEnd,:),1);
        P_new(p).likelihood(1,1)=compute_likelihood(particle, nFirings, blockSize, dt);
        P_new(p).g(1,1) = P_old(p).w((1),1)*P_new(p).likelihood(1,1);
        norm_sum = norm_sum + P_new(p).g(1,1);
    end
    % normalize
    prob_vector = zeros(1, nParticles);
    for p = 1:nParticles
        P_new(p).g(1,1) = P_new(p).g(1,1)/norm_sum;
        prob_vector(1,p) = P_new(p).g(1,1);
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % Step 3: resample % % % % % % % % % % % % %
    resampleIndicies = randsample((1:nParticles), nParticles, 'true', prob_vector);
    for p = 1:nParticles
        P_new(p).theta_p(1,:)=temp.theta_p(resampleIndicies(p),:);
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % Ste 4: draw random samples  % % % % % % %
    for p = 1:nParticles
        flag = true;
        while flag
            v = P_new(p).theta_p(1,:)+ sd.theta_p2*randn(1, nNeurons);
            if v > 0 & v < pi
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
        P_new(p).resampled_likelihood(1,1) = compute_likelihood(particle, nFirings, blockSize, dt);
        P_new(p).w(1,1) = P_new(p).resampled_likelihood(1,1)/P_new(p).likelihood(1,1);
        norm_sum = norm_sum + P_new(p).w(1,1);
    end
    % normalize
    for p = 1:nParticles
        P_new(p).w(1,1) = P_new(p).w(1,1)/norm_sum;
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % Step 6: Compute weighted estimate % % % % % % % % % % %
    for p = 1:nParticles
        x_estimate.theta_p(k,:) = P_new(p).theta_p(1,:)*P_new(p).w(1,1) + x_estimate.theta_p(k,:);
    end
    P_old = P_new;
end

tElapsed = toc;
disp(['Elapsed time = ', num2str(tElapsed/60), 'minutes']);
t2 = 0:dt*blockSize:tmax-dt;
figure, plot(t2, x_estimate.theta_p(:,1), t , theta_p(:,1), t, s, 'LineWidth',2),...
   ylim([0 pi]), grid on;
save results