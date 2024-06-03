% given XX, create normalizex and unnormalizeX min and max values

numSims = size(XX,1);
numx = 7;

xmin_all = zeros(numSims, numx);

for k=1:numSims
    xmin_local = XX{k,2};

    xmin_all(k,:) = xmin_local;
end

normalize_minmax = zeros(2,numx);
normalize_minmax(1,:) = min(xmin_all);
normalize_minmax(2,:) = max(xmin_all);

save("normalize_minmax.mat", "normalize_minmax")