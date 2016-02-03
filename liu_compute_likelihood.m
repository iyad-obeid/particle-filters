function g = liu_compute_likelihood(gInit,sMuMat, local_NMat, dt,nNeurons)

% lambda is [nParticles , blocksz] although all values on each row are the
% same. This is because neither posn nor x change over the duration of the
% block - there may be a more efficient coding implementation of this -
% might want to revisit later (iObeid 5/21/09)

g = gInit;
g = reshape(g,[],1);

for j = 1:nNeurons
    lambda = exp(sMuMat.alpha{j}-((sMuMat.posn - sMuMat.mu{j}).^2 ./ (2*sMuMat.sigma{j}.^2))) * dt;


    lambda_1 = lambda;
    lambda_1(local_NMat{j}~=0) = 0;
    g = g .* prod(exp(-lambda_1),2);

    lambda_1 = lambda;
    lambda_1(local_NMat{j}~=1) = 0;
    lambda_2 = lambda;
    lambda_2(local_NMat{j}~=1) = 1;
    g = g .* prod(exp(-lambda_1) .* lambda_2,2);
end
