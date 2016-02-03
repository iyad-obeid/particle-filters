function [sMu,prior] = liu_vpf_step1(x)
global nNeurons nParticles

A = 1;
B = 20;

sim.mu = normrnd(0,A,nParticles,nNeurons);
sim.posn = normrnd(0,B,nParticles,1);

sMu.mu = x.mu + sim.mu;
sMu.posn = x.posn + sim.posn;

prior.mu   = normpdf(sim.mu,0,A);
prior.posn = normpdf(sim.posn,0,B);

prior = prod(prior.mu,2) .* prior.posn;
