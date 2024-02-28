% given XX, create normalizex and unnormalizeX min and max values

numSims = size(XX,1);

xmin_all = zeros(numSims, 6);

for k=1:numSims
    xmin_local = XX{k,2};

    xmin_all(k,:) = xmin_local;
end

normalize_minmax = zeros(2,6);
normalize_minmax(1,:) = min(xmin_all);
normalize_minmax(2,:) = max(xmin_all);

save("normalize_minmax.mat", "normalize_minmax")