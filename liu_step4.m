function [xNew] = liu_step4(sThetaNew,x)
global nNeurons blocksz nParticles pNames

stdev_stage2.mu = 1;
stdev_stage2.posn = 20;

for i = pNames
    ii = i{1};
    if strcmp(ii,'mu')
        for j = 1:nNeurons
           xNew.mu(:,j) = normrnd(x.mu(:,j)+sThetaNew.mu(:,j) , stdev_stage2.mu);
        end
    else
        eval(sprintf('xNew.%s = normrnd(x.%s+sThetaNew.%s,stdev_stage2.%s);',ii,ii,ii,ii));
    end
end

