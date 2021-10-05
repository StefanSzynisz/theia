function [matrixRandVectors] = correlRandVectors(distributionType, B, numOfVars)
% correlRandVectors = correlRandomVectors(distributionType, B, numOfVar)
%   generates vectors of correlated variables, given
%   distributionType = 'normal' or 'uniform'
%   B = matrix with the correlation structure
%   numOfVars = number of stacked variables in B, eg. 2, 3, 4, etc.

    n_elements = size(B,1);
    if strcmp(distributionType,'normal')
        %------- Normal distribution ----------------------------------------
        phi_uncorrel = randn(n_elements,1); % Generate vector of uncorelated normal randN numbers
        phi_correl = B *phi_uncorrel;  % normal correlated verctor
        
    elseif strcmp(distributionType,'uniform')
        %--- Uniformely distributed random number in the interval(0,1)-----
        phi_uncorrel = rand(n_elements,1); % Generate vector of uncorelated random numbers
        phi_correl = B *phi_uncorrel;
    end
    
    % Chop the long vector phi_correl into variables
    n_rows = n_elements / numOfVars;
    n_cols = numOfVars;
    matrixRandVectors = reshape(phi_correl,[n_rows,n_cols]);
end
