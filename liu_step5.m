function wt = liu_step5(N,blockIndex,x,xNew,sMu,dt)
global nNeurons blocksz nParticles pNames

subBlock_start = (blockIndex-1)*blocksz+1;
subBlock_end = (blockIndex)*blocksz;

[wDim1 , wDim2] = size(sMu.posn); % fix this so it doesn't depend on the name "mu"
% blergh

for j = 1:nNeurons
    local_N       = N(j,subBlock_start : subBlock_end);
    local_N       = reshape(local_N,1,[]);
    local_NMat{j}    = repmat(local_N,nParticles,1);
end

for i = pNames

    ii = i{1};
    
    
    if strcmp(ii,'mu')
        for j = 1:nNeurons
            xNewMat.mu{j} = repmat(xNew.mu(:,j),1,length(local_N));
            sMuMat.mu{j}  = repmat(sMu.mu(:,j),1,length(local_N));
        end
    else
        eval(sprintf('xNewMat.%s = repmat(xNew.%s,1,length(local_N));',ii,ii));
        eval(sprintf('sMuMat.%s = repmat(sMu.%s,1,length(local_N));',ii,ii));
    end

end

for j = 1:nNeurons
    sMuMat.alpha{j} = x.alpha(j) * ones(nParticles,blocksz);
    sMuMat.sigma{j} = x.sigma(j) * ones(nParticles,blocksz);
    xNewMat.alpha{j} = x.alpha(j) * ones(nParticles,blocksz);
    xNewMat.sigma{j} = x.sigma(j) * ones(nParticles,blocksz);
end

numerator   = liu_compute_likelihood(ones(nParticles,1),xNewMat,local_NMat,dt,nNeurons);
denominator = liu_compute_likelihood(ones(nParticles,1),sMuMat ,local_NMat,dt,nNeurons);

wt = numerator./denominator;
wt(denominator == 0 ) = 0; % never divide by zero!

% rescale the quotient weights:

if(sum(wt)>0)
    wt=wt/sum(wt);
else
    wt=ones(nParticles,1)/nParticles; 
    disp('Error in liu_step5.m -> sum of weights is zero');
end

wt = reshape(wt,wDim1,wDim2);


    