function [VAR_1_Rand, VAR_2_Rand] = spatially_correlated_two_variables(VAR_1_MEAN, VAR_1_COV, VAR_2_MEAN, VAR_2_COV, gamma, beta, elem_centroids)

% Statistical input:
% -------------------
% Use consistent units in computations:
% MASS	LENGTH	TIME   FORCE   STRESS	ENERGY	 		
%    g	    mm	  ms	   N	  MPa	  N-mm	
%
% steel      steel           
% DENSITY    YOUNG's    GRAVITY   56.33 km/h = 35 mph
% 7.83e-03  2.07e+05  9.806e-03   15.65 mm/ms
%
n_elements = size(elem_centroids,1);                                        % Number of elements
nVars = 2;                                                                  %two spatial variables are generated
precision_flag = 'single';                                                  %'double' or 'single
plotting_flag = 'yes' ;                                                     % 'yes' or 'no'. 
%plotting works only for prescribed rectangular mesh
%===============

%% Compute the correlation matrices
% based on the distance between the elements and spatial and cross-correlations

% 1) SPATIAL correlation materix:
% based on distance between each element and spatial correlation
% parameter, gamma:
Kspatial = correlMatrix(elem_centroids, gamma, 1.0); 

% 2) Cross-correlation matrix between the variables:
%
Kcross = correlMatrix(elem_centroids, gamma, beta); 

% 3) Assemble the global correlation matrix:
%
Kglobal = [ Kspatial Kcross; Kcross Kspatial];
clear Kspatial Kcross;  % To release memory

%% Create a vector of correlated, random variables
% 4) Eigenvalue diagonalization:
%
if strcmp(precision_flag,'single')
    [Veig,D,~] = eig(Kglobal);                                              % works on both single and double precision.
elseif strcmp(precision_flag,'double')
    [Veig,D,~] = eig(Kglobal);                                              % works on both single and double precision.
    %[Veig,D] = eigs(Kglobal,n_elements);                                   %numerical eigenanalysis, works only on double precision
end
clear Kglobal;                                                              % To release memory

% 5) Compute matrix with the correlation structure:
%
B = Veig * sqrt( abs(D) );
clear Veig D;                                                               % again clear variables to release memory
disp('Eigenvector diagonalization was succesful!'); 

%% Generate random variable:
% 6) Use matrix B with the correlation structure and uncorrelated random
% vectors to generate correlated random vectors.
distributionType = 'normal';                                                % mean=0, std=1, 
% *************************
% or 'uniform' with limits=(-1,1)
matrixRandVectors = correlRandVectors(distributionType, B, nVars);
timestamp = datestr(now,'HHMMSS_FFF');                                      % get_miliseconds to create unique realizations

%% Scale normal random vectors with appropriate mean and std or limits:
VAR_1_Rand = scaleRandVect(matrixRandVectors(:,1), distributionType, VAR_1_MEAN, VAR_1_COV);
min_VAR_1 = min(VAR_1_Rand); max_VAR_1 = max(VAR_1_Rand);

VAR_2_Rand = scaleRandVect(matrixRandVectors(:,2), distributionType, VAR_2_MEAN, VAR_2_COV);
min_VAR_2 = min(VAR_2_Rand); max_VAR_2 = max(VAR_2_Rand);

%% Plot random variables (optional):
if strcmp(plotting_flag,'yes')

    elem_size = norm(elem_centroids(2,2:4)-elem_centroids(1,2:4));
    coordX = elem_centroids(:,2);
    coordY = elem_centroids(:,3);
    
    % ****************************************
    resolution = 1.0 * elem_size;                                           %mm
    saveDataPath = strcat('\data\',timestamp,'\');
    mkdir(saveDataPath);
    %-------------------------------------

    % Young modulus:
    plotTitle = strcat('VAR_1--correl-',' ',num2str(gamma),'--beta-',num2str(beta));
    figure(1); left_pos = 3; bott_pos = 12;                                 %cm
    % Add axis option (manual tuning):
    surfPlot(coordX,coordY,VAR_1_Rand,plotTitle,'', saveDataPath,resolution,...
        left_pos,bott_pos, min_VAR_1, max_VAR_1, 'colorbar_on');

    % Poisson ratio:
    plotTitle = strcat('VAR_2--correl-',' ',num2str(gamma),'mm','--beta-',num2str(beta));
    figure(2); left_pos = 13; bott_pos = 12;                                %cm
    surfPlot(coordX,coordY,VAR_2_Rand,plotTitle,'', ...
        saveDataPath,resolution,left_pos,bott_pos, min_VAR_2, max_VAR_2, 'colorbar_on');
end

%% FUNCTIONS:
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

function scaledRandVect = scaleRandVect(normalRandVect, distributionType, var1, var2)
    
    phi_correl = normalRandVect;    % correlated normal random vector
    n_elements = size(normalRandVect,1);

    if strcmp(distributionType,'normal')
        %------- Normal distribution ----------------------------------------
        variableMEAN = var1;
        variableSTD = var2*variableMEAN;
        
        scaledRandVect = phi_correl * variableSTD + ones(n_elements,1)*variableMEAN;
        
    elseif strcmp(distributionType,'uniform')
        %--- Uniformely distributed random number in the interval(0,1)-----
        variableLowerLIMIT = var1;
        variableUpperLIMIT = var2;
        
        scaledRandVect = phi_correl * (variableUpperLIMIT - variableLowerLIMIT) + variableLowerLIMIT;
    end

end

function matrixK = correlMatrix(elem_centroids, gamma, beta)
% computes correlation matrix
% gamma = spatial correlation scale/length
% beta = cross-correlation
%
    n_elements = size(elem_centroids,1);
    matrixK = zeros(n_elements,n_elements);  % correlation matrix
    for i=1:n_elements
        for j=1:n_elements
            x_i = elem_centroids(i,2);
            y_i = elem_centroids(i,3);
            z_i = elem_centroids(i,4);
        
            x_j = elem_centroids(j,2);
            y_j = elem_centroids(j,3);
            z_j = elem_centroids(j,4);
        
            r_ij = sqrt( (x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2  );
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

end
