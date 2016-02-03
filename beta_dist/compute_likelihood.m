function f = compute_likelihood(particle, N, blockSize, k, dt)
mu = particle.mu(k,:)/pi;
rho = particle.rho(k,:);
sigma_squared = particle.sigma_squared(k,:);

alpha = (mu.^2 - mu.^3 - mu.*sigma_squared)./sigma_squared;
beta = (mu - 2*mu.^2 + mu.^3 - sigma_squared + mu.*sigma_squared)./sigma_squared;

lambda = rho.*betapdf(particle.s(k,1)/pi, alpha, beta);

neuronal_likelihoods = ((lambda*dt).^N).*exp(-lambda*dt*blockSize);

f = prod(neuronal_likelihoods);