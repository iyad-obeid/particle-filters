function [sTheta,sMu,x,rpfx,m] = ...
    liu_step2b(g,sTheta,sMu,x,rpfx,m)

global nNeurons blocksz nParticles pNames

% This function is the resampling portion of step 2

resample_indices = randsample(1:nParticles,nParticles,true,g);
resample_indices = reshape(resample_indices,nParticles,1);

for j = {'sTheta','sMu','x','rpfx','m'}
    jj = j{1};

    for i = pNames
        ii = i{1};

        if (strcmp(jj,'rpfx') && strcmp(ii,'mu'))
            for k = 1:nNeurons
               rpfx.mu{k} = rpfx.mu{k}(resample_indices,:);
            end
        else % normal case
            eval(sprintf('%s.%s = %s.%s (resample_indices,:);',jj,ii,jj,ii));
        end
    end
end

