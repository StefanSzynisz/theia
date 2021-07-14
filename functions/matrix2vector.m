function vector = matrix2vector(matrix,n_rows,n_cols)
% matrix2vector(matrix,n_rows,n_cols)
% convertes matrix into specially ordered vector
    vector = zeros(n_rows*n_cols,1);
    for i = 1:1:n_rows
        matr_row = n_rows+1-i;
        vect_indx = 1 + (i-1)*n_cols;
        vect_indx_end = vect_indx + (n_cols -1);
    
        vector(vect_indx:vect_indx_end) = matrix(matr_row,:);
    end
end