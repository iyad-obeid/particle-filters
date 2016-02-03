function [sThetaNew] = liu_step3(sTheta,m)
global nNeurons blocksz nParticles pNames hh

% sThetaNew = normrnd(m,sqrt(hh)*std(sTheta));

for i = pNames
   ii = i{1};

   % blergh
   if strcmp(ii,'mu')
       for j = 1:nNeurons
           sThetaNew.mu(:,j) = normrnd(m.mu(:,j) , sqrt(hh)*std(sTheta.mu(:,j)));
       end
   else
          eval(sprintf('sThetaNew.%s = normrnd(m.%s,sqrt(hh)*std(sTheta.%s));',ii,ii,ii));
   end
   
end
