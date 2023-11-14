A = 1:5;
B = [5:5:15];
C = 1:3;
num_combs = length(A) * length(B)

[X, Y] = ndgrid(A, B)
X1 = reshape(X, 1, [])
Y1 = reshape(Y, 1, [])

[X1; Y1]

expand_grid_test(A, B, C)

grid = expand_grid(A, B)

get_row_index(grid, 2, 10)

% find(ismember(grid, cell2mat(q), 'rows'))