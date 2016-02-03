function [rpfx] = liu_assign_next_particle(rpfx,xNew,blockIndex)
global nNeurons blocksz nParticles pNames

for i = pNames
    ii = i{1};
    
    if strcmp(ii,'mu')
       for j = 1:nNeurons
           rpfx.mu{j}(:,blockIndex+1)     = xNew.mu(:,j);
       end
    else
        eval(sprintf('rpfx.%s(:,blockIndex+1)     = xNew.%s;',ii,ii));
    end
end
