% Verification problem generator 
%--------------------------------------------------------------------------
% Author: Stefan Szyniszewski
% Date:   06/07/2021
% Description: Script to generate stochastic fields of elastic constants
% and compute corresponding stress and strain fields under given
% boundary conditinos and stress field.
%--------------------------------------------------------------------------
% MAIN
%--------------------------------------------------------------------------
% Input(s):
%   - 
%--------------------------------------------------------------------------
% Ouput(s):
%   - .\results\rnd_mat_prop.mat - random elastic constants 
%  (spatially and cross-correlated)
%   - .\results\epsilon-xx.png, etc - images with results
%--------------------------------------------------------------------------
% See also:
% two_spatially_correl_stochast_vars - function for generation of two 
% correcalted stochastic variables
% LEfe - finite element solver
%--------------------------------------------------------------------------
clearvars; close all; clc; clear; tic; 
[thisPath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(thisPath);  addpath('../functions/');                                     % add path to functions folder                                                             % change directory to current path

%% Specify key parameters for FE solver function:
% Specify dimmensions:
% Use consistent units in computations:                                     % For example, STEEL:       
% MASS	LENGTH	TIME   FORCE   STRESS	ENERGY	                            % DENSITY    YOUNG's    GRAVITY   56.33 km/h = 35 mph		
%    g	    mm	  ms	   N	  MPa	  N-mm	                            % 7.83e-03  2.07e+05  9.806e-03   15.65 mm/ms                            
%
%      Y
%      ^
%      |     pressure load
%      |     -------------
%      | || || || || || || || ||
%      | \/ \/ \/ \/ \/ \/ \/ \/        |
%      |------------------------.       -
%      | 26 | 27 | 28 | 29 | 30 |       |
%      |------------------------|       |
%      | 21 | 22 | 23 | 24 | 25 |       |
%      |------------------------|       |
%      | 16 | 17 | 18 | 19 | 20 |      Ly   (ny)
%      |------------------------|       |
%      | 11 | 12 | 13 | 14 | 15 |       |
%      |------------------------|       |
%      |  6 |  7 |  8 |  9 | 10 |       |
%      |------------------------|       |
%      |  1 |  2 |  3 |  4 |  5 |       |
%      -----------------------------------------> X
%      ^    ^    ^    ^    ^    ^
%       simply supported boundary
%                 
%     -/-------- Lx -----------/--
%    (nx - number of elements in X-direction)
%

% Dimmnesions:
Lx = 100;                                                                   %mm
Ly = Lx;                                                                    %mm, same dimmensions
nx = 50;
ny = nx;                                                                    %same number of subdivsions

% Based on the regular geometry compute the coordinates of element centroids
% such that we can compute the distance between the elements for
% correlation computations 
[coordX, coordY] = elemCoordVector(Lx,Ly,nx,ny);

%% Physical problem description:
% *****************************
pressure = -1.0;                                                            % pressure in [MPa]
bc_type = 'pressure';                                                       % type of boundaries at the bottom: 'fixed' or 'roller' or 'pressure' 

%% Folder for saving the results:
timestamp = datestr(now,'HHMMSS_FFF');                                      % get_miliseconds to create unique realizations
saveDataFolder = strcat(thisPath,'\data\',timestamp,'\');                   % data folder
if ~exist(saveDataFolder, 'dir')
    mkdir(saveDataFolder)                                                   % check if exists
end

%% Statistical parameters:
% -----------------------
E_mean = 1000;                                                              % MPa, 3D printed composite
E_cov = 0.2;                                                                % coefficient of variation
v_mean = 0.4;
v_cov = 0.1;                                                                % coefficient of variation

distributionType = 'normal';                                                % mean=0, std=1, 
% or 'uniform' with limits=(-1,1)
%==========================

%% Compute the stochastic field:
% gamma = 70;                                                                 % mm, spatial correlation length
% beta = 0.5;                                                                 % varying from 0.0 to 1.0 is the cross correlation between the two variables
% 
% matrixRandVectors = generate_spatially_correl_vars(coordX, coordY, gamma, beta, distributionType);
% 
% % Save the random variables into .mat file;
% rndFileName = [saveDataFolder 'matrixRandVectors_' timestamp '.mat'];
% save(rndFileName,'matrixRandVectors');                                     % save unit correlated random variables
% disp(['Random variables were written to' ' ' rndFileName]);

%% Read stochastic variable from the previous run:
read_timestamp = '154100_122';                                              % resuse existing random field
load_from_file = ['data\' read_timestamp '\matrixRandVectors_' read_timestamp '.mat'];
tmp = load(load_from_file,"matrixRandVectors");                             % structure
matrixRandVectors = tmp.matrixRandVectors;                                  % extract from the structure   

%% Scale normal random vectors with appropriate mean and std or limits:
E_rand = scaleRandVect(matrixRandVectors(:,1), distributionType, E_mean, E_cov);
E_min = min(E_rand); E_max = max(E_rand);

v_rand = scaleRandVect(matrixRandVectors(:,2), distributionType,v_mean, v_cov);
v_min = min(v_rand); v_max = max(v_rand);

%% Plot synthetic stochastic material parameters:
fig_num = 0;                                                                % initialize figure number
% Elastic constants:
fig_num = fig_num +1; %close(fig_num);                                      % figure number
my_subplot(fig_num,{'Young', 'Poisson'},[E_rand,v_rand],nx, ny);
saveas(gcf,[saveDataFolder 'material-benchmark-' bc_type timestamp '.png']);     % save images as png

%% Set FE parameters based on the user input and characteristics of the supplied data:
mesh = setupmesh(Lx,Ly,nx,ny,pressure,E_rand,v_rand,bc_type);               % initialize mesh

%% Solve FE problem:
[sig,eps,~] = LEfe(mesh);                                                  % finite element solver

% eps_xx,eps_yy,eps_yy,sig_xx,sig_yy,sig_xy,E,v
% --- strain ---, --- stress ---,--- mat ---
results = [eps,sig,E_rand,v_rand]; 
result_file_path = [saveDataFolder 'benchmark-' bc_type '_' timestamp '.xlsx'];
writematrix(results, result_file_path);                                     % write results to Excel file    

%% Plot results:
% Strains:
plot_titles = {'epsilon-xx - benchamrk', 'epsilon-yy - benchmark', 'epsilon-xy - benchmark'};
fig_num = fig_num +1; %close(fig_num);                                      % figure number
my_subplot(fig_num,plot_titles,eps,nx,ny);
saveas(gcf,[saveDataFolder 'strains-benchmark-' bc_type '_' timestamp '.png']);     % save images as png

% Stress:
fig_num = fig_num +1; %close(fig_num);                                      % figure number
my_subplot(fig_num,{'sigma-xx - benchmark', 'sigma-yy - benchmark', 'sigma-xy - benchmark'},sig,nx,ny);
saveas(gcf,[saveDataFolder 'stresses-benchmark-' bc_type '_' timestamp '.png']);    % save images as png


%% FUNCTIONS:
function [coordX, coordY] = elemCoordVector(Lx,Ly,nx,ny)
    n_elements = nx * ny;
    deltaX = Lx / nx;  % size of element in x-direction
    deltaY = Ly / ny;  % size of element in y-direction
    
    coordX = zeros(n_elements,1); % placeholder for the element coordinates
    coordY = zeros(n_elements,1); % 
    
    %index of the element position in the grid
    ind_gridX =1;  % start at the first element in position (1,1)
    ind_gridY =1;
    
    %   ind_gridY
    %      ^
    %      |                             
    %      |-----------------------------.   
    %      | 6,1 | 6,2 | 6,3 | 6,4 | 6,5 | 
    %      |-----------------------------| 
    %      | 5,1 | 5,2 | 5,3 | 5,4 | 5,5 | 
    %      |-----------------------------| 
    %      | 4,1 | 4,2 | 4,3 | 4,4 | 4,5 | 
    %      |-----------------------------|       
    %      | 3,1 | 3,2 | 3,3 | 3,4 | 3,5 |  
    %      |-----------------------------|  
    %      | 2,1 | 2,2 | 2,3 | 2,4 | 2,5 |       
    %      |-----------------------------|       
    %      | 1,1 | 1,2 | 1,3 | 1,4 | 1,5 |       
    %      ------------------------------------> ind_gridX
    %                            (nx * ny)
    %                  

    for elem_id=1:n_elements  % march element by element:
        % as we move from element to the next element
        coordX(elem_id,1) = deltaX/2 + (ind_gridX -1)* deltaX; % increment x-coordinate
        
        % as we move from the row to the next row of elements
        coordY(elem_id,:) = deltaY/2 + (ind_gridY -1)* deltaY; % increment y-coordinate
    
        ind_gridX = ind_gridX+1;  % increment index x as we progress
        if ind_gridX > nx  % if we exceed the number of subdivisions in x-direction
            % then reset the ind_gridX=1 and increment ind_gridY by +1
            ind_gridX = 1;
            ind_gridY = ind_gridY +1;
        end
    end
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
