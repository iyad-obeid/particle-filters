function g = liu_step2a(x,N,w,blockIndex,sMu,dt)
global nNeurons blocksz nParticles pNames

subBlock_start = (blockIndex-1)*blocksz+1;
subBlock_end = (blockIndex)*blocksz;

[wDim1 , wDim2] = size(w);

% local_posn    = posn(subBlock_start : subBlock_end);

% local_posn    = reshape(local_posn,1,[]);
% local_posnMat = repmat(local_posn,nParticles,1);

for j = 1:nNeurons
    local_N       = N(j,subBlock_start : subBlock_end);
    local_N       = reshape(local_N,1,[]);
    local_NMat{j}    = repmat(local_N,nParticles,1);
end

for i = pNames

    if strcmp(i,'mu')
        for j = 1:nNeurons
%             eval(sprintf('sMu.%s(:,j) = reshape(sMu.%s,[],1);',i{1},i{1}));
            eval(sprintf('sMuMat.%s{j} = repmat(sMu.%s(:,j),1,blocksz);',i{1},i{1}));
        end
    else
        eval(sprintf('sMu.%s = reshape(sMu.%s,[],1);',i{1},i{1}));
        eval(sprintf('sMuMat.%s = repmat(sMu.%s,1,blocksz);',i{1},i{1}));
    end

end

for j = 1:nNeurons
    sMuMat.alpha{j} = x.alpha(j) * ones(nParticles,blocksz);
    sMuMat.sigma{j} = x.sigma(j) * ones(nParticles,blocksz);
end

g = liu_compute_likelihood(w,sMuMat,local_NMat,dt,nNeurons);

% rescale weights
if(sum(g)>0)
    g=g/sum(g);
else
    g=ones(nParticles,1)/nParticles;
    disp('Error in liu_step2.m -> sum of weights is zero');
end

g = reshape(g,wDim1,wDim2);