function matrixK = correlMatrix(coordX, coordY, gamma, beta)
% computes correlation matrix
% gamma = spatial correlation scale/length
% beta = cross-correlation
%
    n_elements = size(coordX,1);
    matrixK = zeros(n_elements,n_elements);  % correlation matrix
    for i=1:n_elements
        for j=1:n_elements
            x_i = coordX(i,1);
            y_i = coordY(i,1);
        
            x_j = coordX(j,1);
            y_j = coordY(j,1);
        
            r_ij = sqrt( (x_i - x_j)^2 + (y_i - y_j)^2  );
            k_ij = beta *exp( -r_ij^2 / gamma^2 );
            % round small numbers to zero
            tolerance = 1e-4;
            if (k_ij<tolerance)
                k_ij=0;
            end
            matrixK(i,j) = k_ij;
            
        end
    end
end