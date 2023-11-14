function output_grid = expand_grid(varargin)
    % EXPAND_GRID Creates a matrix with all combinations of the vectors
    % supplied as arguments.

    % Creates our output, but stored as a separate tensor for each
    % parameter
    [output_mat{1:nargin}] = ndgrid(varargin{:});

    % Convert from a list of tensors to a tall matrix
    for i=1:length(output_mat)
        this_mat = output_mat{i};
        % Convert tensor to a long vector
        this_vec = reshape(this_mat, 1, []);
        % Store vector as a column of the output matrix
        output_grid(:,i) = this_vec; 
        % Matlab wants me to pre-alloate the output matrix. I don't want
        % to unless the function runs really slowly. I don't expect this to
        % be the case, but if it is I'm happy to adjust.
    end   
end