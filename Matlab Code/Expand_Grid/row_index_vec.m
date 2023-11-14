function ind = row_index_vec(comb_grid, par_vec)
% ROW_INDEX_VEC Gets the row index in comb_grid of the row specified by
% par_vec (order of entries in par_vec matters!).
    ind = find(ismember(comb_grid, par_vec, 'rows'));
end
    