function f = compute_likelihood(particle, N, blockSize, dt, nNeurons)

lambda = particle.beta(1,:) + particle.gamma(1,:).*cos(2*(particle.s(1,1) - particle.theta_p(1,:)));

for n = 1:nNeurons
    if particle.theta_p(1,n) > pi/2 && ...
            particle.theta_p(1,n)-pi/2 < particle.s(1,1)
        % lambda(1,n)=1e-10;
        N(n)=0;
    end
    if particle.theta_p(1,n) < pi/2 && ...
            particle.theta_p(1,n)+pi/2 > particle.s(1,1)
        % lambda(1,n)=1e-10;
         N(n) = 0;
    end
end
neuronal_likelihoods = ((lambda*dt).^N).*exp(-lambda*dt*blockSize);
f = prod(neuronal_likelihoods);
