%% Display good_indices parameters

numGood = length(good_indices);

numInput = 18;

good_inputs = zeros(numGood,numInput); % 18 input

for k =1:length(good_indices)
    current_index = good_indices(k);
    paramk = XX{k,1};
    good_inputs(k,:) = paramk;
end