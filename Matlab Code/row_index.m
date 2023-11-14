function ind = row_index(comb_grid, varargin)
    % ROW_INDEX Gets the row index in comb_grid of the row corresponding to
    % the values supplied as arguments after comb_grid.
    % Note: The values supplied after comb_grid must collectively specify
    % exactly one row of comb_grid (in the correct order).
    ind = find(ismember(comb_grid, cell2mat(varargin), 'rows'));
end
    