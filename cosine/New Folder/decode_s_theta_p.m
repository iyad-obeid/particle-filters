clc
clear all
close all
load input_data

nParticles = 100;
blockSize  = 10;
nIterations= floor(length(t)/blockSize);

sd.theta_p1 = 0.0001;
sd.theta_p2 = 0.0001;

sd.s1 = 0.05;
sd.s2 = 0.05;

gamma_const = gamma_start;
beta_const = beta_start;

P_new(1:nParticles) = struct('s',0,...
    'gamma', gamma_const.*ones(1, nNeurons),...
    'theta_p',zeros(1, nNeurons),...
    'beta', beta_const.*ones(1, nNeurons),...
    'g',zeros(1,1),...
    'likelihood',zeros(1,1),...
    'resampled_likelihood',zeros(1,1),...
    'w',ones(1,1));
P_old = P_new;

s_ds = downsample(s, blockSize);
x_estimate.theta_p= zeros(nIterations, nNeurons);
x_estimate.s= zeros(nIterations, 1);
% initialize first iteration
for p =1:nParticles
    P_old(p).theta_p(1,:)= theta_p(1,:);
    P_old(p).s(1,1)=s_ds(1,1);
end
x_estimate.theta_p(1, :) = theta_p_start;
x_estimate.s(1,1)= s_ds(1,1);


for k = 2:nIterations
    disp(['block number ', num2str(k), ' of ', num2str(nIterations)]);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % Step1: draw random samples  % % % % % % %

    for p = 1:nParticles
        flag = true;
        while flag
            v = P_old(p).s(1,1)+ sd.s1*randn;
            if v > 0 && v < pi
                P_new(p).s(1,1) = v;
                flag = false;
            end
        end
        flag = true;
        while flag
            v = P_old(p).theta_p(1,:)+ sd.theta_p1*randn(1, nNeurons);
            if v > 0 & v < pi
                P_new(p).theta_p(1,:) = v;
                flag = false;
            end
        end
        temp.theta_p(p,:) = P_new(p).theta_p(1,:);
        temp.s(p,:) = P_new(p).s(1,1);
    end

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % Step2: compute likelihoods % % % % % % % %
    norm_sum = 0;
    for p = 1:nParticles
        particle = P_new(p);
        blockStart = (k-1)*blockSize + 1;
        blockEnd   = blockStart + blockSize - 1;
        nFirings = sum(N(blockStart:blockEnd,:),1);
        P_new(p).likelihood(1,1)=compute_likelihood(particle, nFirings, blockSize, dt, nNeurons);
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
        P_new(p).s=temp.s(resampleIndicies(p),1);
    end

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % Ste 4: draw random samples  % % % % % % %
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
        P_new(p).resampled_likelihood(1,1) = compute_likelihood(particle, nFirings, blockSize, dt, nNeurons);
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
        x_estimate.s(k,1) = P_new(p).s(1,1)*P_new(p).w(1,1) + x_estimate.s(k,1);
    end
    P_old = P_new;
end

t2 = 0:dt*blockSize:tmax-dt;
figure, plot(t2, x_estimate.s, t2, s_ds);

save results
