function [sMu,sTheta,m]=liu_step1(x)
global nNeurons blocksz nParticles pNames hh

a = sqrt(1-hh^2);

for i = pNames

    sTheta_local      = zeros(nParticles,1);
    eval(sprintf('x_local = x.%s;',i{1}));

    if strcmp(i,'mu')
        sTheta_local = normrnd(0,1,nParticles,nNeurons);
    elseif strcmp(i,'posn');
        sTheta_local = normrnd(0,20,nParticles,1); % 60 is a guess based on posn slope
    end

    m_local     = a*sTheta_local + repmat( (1-a)*mean(sTheta_local,1), nParticles , 1);                 % eqn 10.3.9
    theta_local = m_local        + hh*repmat(std(sTheta_local,1), nParticles, 1)...
        .* randn(size(sTheta_local)); % eqn 10.3.8
    sMu_local   = x_local        + theta_local;

    eval(sprintf('sTheta.%s = sTheta_local;',i{1}));
    eval(sprintf('m.%s      = m_local;',i{1}));
    eval(sprintf('sMu.%s    = sMu_local;',i{1}));
    clear sTheta_local m_local sMu_local
end

