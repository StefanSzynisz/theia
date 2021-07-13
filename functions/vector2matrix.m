function matrix = vector2matrix(vector,n_rows,n_cols)
% vector2matrix(vector,n_rows,n_cols)
% converts vector to matrix:
    matrix = zeros(n_rows, n_cols); % Pre-allocate the matrix
    for i = 1:1:n_rows
        matr_row = n_rows+1-i;
        vect_indx = 1 + (i-1)*n_cols;
        vect_indx_end = vect_indx + (n_cols -1);
        matrix(matr_row,:) = vector(vect_indx:vect_indx_end);
    end
end