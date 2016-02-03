function [xEstimate] = liu_assign_next_estimate(xEstimate,w,rpfx,blockIndex)
global nNeurons blocksz nParticles pNames

for i = pNames
    ii = i{1};
    if strcmp(ii,'mu')
        for j = 1:nNeurons
            xEstimate.mu(blockIndex+1,j) = sum(w .* rpfx.mu{j}(:,blockIndex));
        end
    else
        eval(sprintf('xEstimate.%s(blockIndex+1) = sum(w .* rpfx.%s(:,blockIndex));',ii,ii));
    end
end
