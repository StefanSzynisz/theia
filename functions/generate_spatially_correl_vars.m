function [matrixRandVectors] = generate_spatially_correl_vars(coordX, coordY, gamma, beta, distributionType)
%Two dimensional generator of spatially correlated, stochastic variables
%--------------------------------------------------------------------------
% Author: Stefan Szyniszewski
% Date:   05/10/2021
% Description:
% Function to generate a spatial field of two spatially correlatied, 
% and cross correlated variables.
%--------------------------------------------------------------------------
% [a b] = two_spatially_correl_stochast_vars(x,y,z)
%--------------------------------------------------------------------------
% Input(s):
% xxx - yyy
% xxx - yyy
%--------------------------------------------------------------------------
% Ouput(s);
% zzz  - xxx
% zzz - xxx
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------
nVars = 2;                                                                  % it is a field of two variables

% 1) using centroid coordinates of each element, compute the 
% SPATIAL correlation matrix based on distance between each element and 
% spatial correlation parameter, gamma:
Kspatial = correlMatrix(coordX, coordY, gamma, 1.0); 
disp('Spatial correlation created.'); 

% 2) Compute cross correlation matrix between the two variables:
Kcross = correlMatrix(coordX, coordY, gamma, beta);
disp('Cross-correlation matrix created.')

% 3) Assemble the global correlation matrix:
Kglobal = [Kspatial Kcross; Kcross Kspatial];
clear Kcross Kspatial;

% 4) Eigenvalue diagonalization:
[Veig,D,~] = eig(Kglobal);
clear Kglobal;                                                              % To release memory
disp('Eigenvector diagonalization was succesful!'); 

% 5) Compute matrix with the correlation structure:
B = Veig * sqrt( abs(D) );
clear Veig D;                                                               % again clear variables to release memory
disp('Matrix with correlation structure created.')

%% Generate random variable:
% 6) Use matrix B with the correlation structure and uncorrelated random
% vectors to generate correlated random vectors.
matrixRandVectors = correlRandVectors(distributionType, B, nVars);
end

