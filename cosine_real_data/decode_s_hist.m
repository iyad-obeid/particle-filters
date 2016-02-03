clc
clear all
%close all
load kinematics
load init_params


start_sample = training_period/dt + 1;
decode_period = 30;

tStart = training_period + dt;
tmax = training_period + decode_period;
t = (0:dt:tmax-dt)';
s = kinematics(1:tmax*(1/dt), 2);
clear kinematics

%figure, plot(t, s)

% N_low = 1000;
% N_high = 15000;
% f = find(N_low < N_total & N_total < N_high);
% neuron_list = neuron_list(f);
% gamma_const = gamma_const(f);
% beta_const = beta_const(f);

neuron_list = 1:158;
% neuron_list = [    02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 ...
%                 21 22 23 24 25    27 28 29 30 31 32    34 35    37 38 39 40 ...
%                 41 42 43 44    47    49 50 51    53 54 55  57 58 59 60 ...
%                 61  63 64    66 67 68 69 70 71    73    75 76 77 78       ...
%                 81       84 85          89       92    94 95   98  99 100 ...
%                        105 107 110 111 112 113 114 115 116 117 118 119 120   ...
%                  122  124  126 127 128 131 132 133 134 135 136 139 140 ...
%                 141 142 143 144 145 149 150 151 152 153 154 155 156 158];  

% neuron_list = neuron_list(r);
s_init = s((start_sample - 1),1);
theta_p_init = pd(neuron_list);
gamma_const  = gamma_const(neuron_list);
beta_const   = beta_const(neuron_list);

nNeurons = length(neuron_list);

% load firing times for neuron 'm'
for m = 1:length(neuron_list)
    neuron_number = neuron_list(m);
    eval(['load neuron', num2str(neuron_number)]);
    eval(['N(:,', num2str(m),') = N', num2str(neuron_number), '(1:length(s));']);
    eval(['clear N', num2str(neuron_number)]);
end


s = s(start_sample:end,1);
N = N(start_sample:end,:);

blockSize = 20;
s_ds = downsample(s, blockSize);

nIterations = floor((tmax-training_period)/(blockSize*dt));
nParticles = 100;
sd.s1 = 0.05;
block number 1498 of 1500
block number 1499 of 1500
sd.s2 = 0.025;
sd.theta_p1 = 0.008;
sd.theta_p2 = 0.001;

P_new(1:nParticles) = struct('s',0,...
    'gamma', gamma_const.*ones(1, nNeurons),...
    'theta_p',zeros(1, nNeurons),...
    'beta', beta_const.*ones(1, nNeurons),...
    'g',zeros(1,1),...
    'likelihood',zeros(1,1),...
    'resampled_likelihood',zeros(1,1),...
    'w',ones(1,1));
P_old = P_new;

% initialize first iteration
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
        temp.s(p,1) = P_new(p).s(1,1);
        
        P_new(p).theta_p(1,:) = theta_p_init;
        %  temp.theta_p(p,:) = P_new(p).theta_p(1,:);
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
        %  P_new(p).theta_p(1,:)=temp.theta_p(resampleIndicies(p),:);
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
        x_estimate.s(k,1) = P_new(p).s(1,1)*P_new(p).w(1,1) + x_estimate.s(k,1);
    end
    P_old = P_new;
end
tElapsed = toc;
disp(['Elapsed time = ', num2str(tElapsed/60), 'minutes']);

t2 = tStart:dt*blockSize:tmax-dt;
L1 = ['# neurons = ', num2str(nNeurons), ...
    ' # particles = ', num2str(nParticles), ...
    ' block size = ', num2str(blockSize)];
L2 = ['\sigma_{s1} = ' ,num2str(sd.s1), '   \sigma_{s2}', num2str(sd.s2)];
L3 = ['Bin size = ', num2str(db)];

figure, plot(t2, x_estimate.s,'b'), grid on, hold on, ...
    plot(t2, s_ds, 'r','LineWidth',2), ...
    xlabel('time (s)'), ylabel('\Theta'), title({L1; L2; L3});


save results

% figure, 
% for k = 1:nNeurons
%     hold on, plot(x_estimate.theta_p(:,k));
% end