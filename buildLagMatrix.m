function m2 = buildLagMatrix(m,lag)

rng = 1:size(m,1) - lag + 1;
nCols = size(m,2);

m2 = zeros(length(rng),lag*nCols + 1);
k = 2;
for i = 1:lag
    for j = 1:nCols
        m2(:,k) = m(rng,j);
        k = k + 1;
    end
    rng = rng + 1;
end

m2(:,1) = 1;
