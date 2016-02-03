clc
clear all
close all
load kinematics
tmax = 20;
t = (0:dt:tmax-dt)';
s = kinematics(2, 1:tmax*(1/dt))' + 2.5;
clear kinematics

% neuron_list = [33, 45, 46, 52, 56, 60, 66, 69, 71, 74, 80, 82, 83, ...
%     86, 87, 88, 89, 91, 92, 93, 95, 97, 98, 99, 101, 102, 103, 106, ...
%     107, 108, 109, 112, 115, 116, 117, 118, 120, 121, 126, 127, 129, ...
%     133, 134, 135, 136, 137, 138, 140, 141, 144, 145, 146, 147, 148, ...
%     149, 150, 151, 152, 153, 155, 156, 157, 158];

neuron_list = 1;

nParticles = 20;
blockSize  = 50;
var_mu1 = 0.0001;
var_mu2 = 0.0001;
rho_const = 15;
sigma_squared_const = 0.05;

nIterations= floor(length(t)/blockSize);
nSimulations = 5;

P(1:nParticles) = struct('s',downsample(s, blockSize),...
    'rho', rho_const*ones(nIterations, 1),...
    'mu',zeros(nIterations, 1),...
    'sigma_squared', sigma_squared_const*ones(nIterations, 1),...
    'g',zeros(nIterations,1),...
    'likelihood',zeros(nIterations,1),...
    'resampled_likelihood',zeros(nIterations,1),...
    'w',ones(nIterations,1));
x_estimate.mu= zeros(nIterations, 1);

for m = 1:length(neuron_list)
    % create directory for neuron
    neuron_number = neuron_list(m);
    mkdir(['neuron', num2str(neuron_number)]);
    % load firing times for neuron 'm'
    eval(['load neuron', num2str(neuron_number)]);
    eval(['N(:,1) = N', num2str(neuron_number), '_resampled(1:length(s));']);
    
    for sim_num = 1:nSimulations
        disp(['simulation #', num2str(sim_num),...
            ' of ', num2str(nSimulations), ' for neuron ', num2str(neuron_number)]);
        P(1:nParticles) = struct('s',downsample(s, blockSize),...
            'rho', rho_const*ones(nIterations, 1),...
            'mu',zeros(nIterations, 1),...
            'sigma_squared', sigma_squared_const*ones(nIterations, 1),...
            'g',zeros(nIterations,1),...
            'likelihood',zeros(nIterations,1),...
            'resampled_likelihood',zeros(nIterations,1),...
            'w',ones(nIterations,1));
        x_estimate.mu= zeros(nIterations, 1);
        
        % initialize first iteration
        for p =1:nParticles
            P(p).mu(1,:)= pi/2;
            x_estimate.mu(1, :) = pi/2;
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
                disp(p);
                disp(['1st stage: ', num2str(P(p).g(k,1))]);
                disp(['2nd stage: ', num2str(P(p).w(k-1,1))]);
                disp(['s: ', num2str(P(p).s(k,1))]);
                disp(['mu: ', num2str(P(p).mu(k,1))]);
                
                particle = P(p);
                blockStart = (k-1)*blockSize + 1;
                blockEnd   = blockStart + blockSize - 1;
                nFirings = sum(N(blockStart:blockEnd,:),1);
                P(p).likelihood(k,1)=compute_likelihood(particle, nFirings, blockSize, k, dt);
                P(p).g(k,1) = P(p).w((k-1),1)*P(p).likelihood(k,1);
                norm_sum = norm_sum + P(p).g(k,1);
                
                
                disp(['likelihood: ', num2str(P(p).likelihood(k))]);

            end
           
            % normalize
            prob_vector = zeros(1, nParticles);
            for p = 1:nParticles
                P(p).g(k,1) = P(p).g(k,1)/norm_sum;
                prob_vector(1,p) = P(p).g(k,1);
            end
            
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % % % % % % % % Step 3: resample % % % % % % % % % % % % %
           % prob_vector(isnan(prob_vector))=0;
            
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
    end
    tElapsed = toc;
    disp(['Elapsed time = ', num2str(tElapsed/60), 'minutes']);
    t2 = 0:dt*blockSize:tmax-dt;
    fname = [cd,'/neuron', num2str(neuron_number), '/simulation', num2str(sim_num)];
    save(fname, 'x_estimate', 't2');
end
save neuron_list neuron_list nSimulations
figure, plot(t2, x_estimate.mu(:,1), t, s, 'LineWidth',2),...
    ylim([0 pi]), grid on;