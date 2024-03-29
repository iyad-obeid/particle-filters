clc
clear all
close all
load kinematics
load initial_estimates

decode_period = (40:dt:(70-dt));
tStart =  decode_period(1);
tEnd = decode_period(end);
t = (tStart:dt:tEnd-dt)';
startSample = (decode_period(1) + dt)/dt;
endSample = (decode_period(end) + dt)/dt;
samples = (startSample:endSample)';
s = kinematics(samples, 2);

clear kinematics blockSize

neuron_list = 1:158;
% neuron_list = [ 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 ...
%     21 22 23    25    27    29 30 31 32    34 35    37    39 40 ...
%     41 42 43 44             49 50 51    53 54 55 56 57 58 59 60 ...
%     61 62 63 64 65 66 67 68 69 70 71    73 74    76 77 78       ...
%           83 84 85                   92 93 94 95 96 97 98 99 100 ...
%     101  105  110 111  113 114 115 116  119 120   ...
%     122 123 124 125 126 128 131 132 133 135 136 137 139 140 ...
%     141 142 143 145       150 151 152 153 154 155 156 ];

% neuron_list = [ 07 08 10 14 15 16 17 18 19 20 ...
%     21 22 23    25 32    34  37   40 ...
%     41 42     49 50 51    53 55 56 57 58 60 ...
%     61 62 63 64 65 66 67 68 69 70 71    73 74    76 77 78       ...
%           83 84 85                   92 93 94 95 97 98 99 100 ...
%     101  105  110 111  113 114 115 116  119 120   ...
%     122 123 126 128 131 132 133 135 136 137 139 140 ...
%     141 142 145       150 151 152 153 155 156 ];

% theta_p_init = pd(neuron_list);
theta_p_init = theta_p_final(neuron_list);
gamma_const  = gamma_final(neuron_list);
beta_const   = beta_final(neuron_list);

nNeurons = length(neuron_list);

% load firing times for neuron 'm'
for m = 1:length(neuron_list)
    neuron_number = neuron_list(m);
    eval(['load neuron', num2str(neuron_number)]);
    eval(['N(:,', num2str(m),') = N', num2str(neuron_number),';']);
    eval(['clear N', num2str(neuron_number)]);
end
% s = s(start_sample:end,1);
N = N(samples,:);

blockSize = 50;
s_ds = downsample(s, blockSize);
nIterations= floor(length(s)/blockSize);
nParticles = 500;
sd.s1 = 0.07;
sd.s2 = 0.02;
sd.theta_p1 = 0.01;
sd.theta_p2 = 0.001;

P_new(1:nParticles) = struct('s',0,...
    'gamma', zeros(1, nNeurons),...
    'theta_p',zeros(1, nNeurons),...
    'beta', zeros(1, nNeurons),...
    'g',zeros(1,1),...
    'likelihood',zeros(1,1),...
    'resampled_likelihood',zeros(1,1),...
    'w',ones(1,1));
P_old = P_new;

% initialize first iteration
s_init = s_ds(1);
for p =1:nParticles
    P_old(p).theta_p(1,:) = theta_p_init;
    P_old(p).s(1,1) = s_init;
end
x_estimate.theta_p= zeros(nIterations, nNeurons);
x_estimate.s= zeros(nIterations, 1);

tic;
for k = 1:nIterations
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
            v = P_old(p).theta_p(1,:)+ sd.theta_p1*randn(1,nNeurons);
            if v > 0 & v < pi
                P_new(p).theta_p(1,:) = v;
                flag = false;
            end
        end
        
        temp.s(p,1) = P_new(p).s(1,1);
        
        % P_new(p).theta_p(1,:) = theta_p_init;
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
        P_new(p).s(1,1)=temp.s(resampleIndicies(p),1);
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
            v = P_old(p).theta_p(1,:)+ sd.theta_p2*randn(1,nNeurons);
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
        x_estimate.s(k,1) = P_new(p).s(1,1)*P_new(p).w(1,1) + x_estimate.s(k,1);
        x_estimate.theta_p(k,:) = P_new(p).theta_p(1,:)*P_new(p).w(1,1) + x_estimate.theta_p(k,:);
        x_estimate.gamma(k,:) = P_new(p).gamma(1,:)*P_new(p).w(1,1) + x_estimate.gamma(k,:);
        x_estimate.beta(k,:) = P_new(p).beta(1,:)*P_new(p).w(1,1) + x_estimate.beta(k,:);
    end
    P_old = P_new;
end
tElapsed = toc;
disp(['Elapsed time = ', num2str(tElapsed/60), 'minutes']);

t2 = tStart:dt*blockSize:tEnd-dt;
L1 = ['# neurons = ', num2str(nNeurons), ...
    ' # particles = ', num2str(nParticles), ...
    ' block size = ', num2str(blockSize)];
L2 = ['\sigma_{s1} = ' ,num2str(sd.s1), '   \sigma_{s2} = ', num2str(sd.s2)];


figure, plot(t2, x_estimate.s,'b'), grid on, hold on, ...
    plot(t2, s_ds, 'r','LineWidth',2), ...
    xlabel('time (s)'), ylabel('\theta (rad)'), ...
    xlim([tStart tEnd]), ylim([0 pi]), ...
    legend('\theta_{est}', '\theta_{true}'), ...
  %  title({L1; L2;});


% save results

figure, ylim([0 pi]), title('\theta est'),
for k = 1:nNeurons
    hold on, plot(x_estimate.theta_p(:,k));
end

figure, title('\gamma est'),
for k = 1:nNeurons
    hold on, plot(x_estimate.gamma(:,k));
end

figure, title('\beta est'),
for k = 1:nNeurons
    hold on, plot(x_estimate.beta(:,k));
end

