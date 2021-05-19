%plane stress elastic property solver 
%--------------------------------------------------------------------------
% Authors: William Coombs, Stefan Szyniszewski
% Last modified:   21/10/2020
% Description:
% Script to determine the properties of a material based on a known
% strain field, epsA, and predicted stress field from linear finite
% element analysis.  The script iterates to find the elastic properties of
% the material and stress field until it converges within a given 
% tolerance or reaches a maximum number of iteraitons (default set at 50). 
%
% The script requires:
% - xlsx file with the known strain information: Strain.xlsx
% - settings.yml with information about the applied load and other
%   parameters to define the physical problem.
%
%
%--------------------------------------------------------------------------
% MAIN
%--------------------------------------------------------------------------
% Input(s):
% 
%--------------------------------------------------------------------------
% Ouput(s);
% 
%--------------------------------------------------------------------------
% See also:
% settings.yml      - settings and description of the problem
% SETUPMESH         - problem set up 
% LEFE              - finite element solver
% UPDATEELASTICPROP - elastic property update
%--------------------------------------------------------------------------
clearvars; close all; clc; clear; tic;
[thisPath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(thisPath); addpath('functions')                                          % add path to functions folder
                                                                            % ************* designates input                                            

%% Read YAML file:
settings = YAML.read('settings.yml');                                       % use function from http://vision.is.tohoku.ac.jp/~kyamagu/software/yaml/
% **********************************
timestamp = datestr(now,'HHMMSS_FFF');                                      % get_miliseconds to create unique realizations

%% Set values of parameters using info from 'settings.yml' file:
load(settings.RndMatPropPath,'YoungRand','PoissonRand');                    % Benchmark data to test the algorithm
epsA = xlsread(settings.StrainDataFilePath);                                % read strain data
sigA = xlsread(settings.StressDataFilePath); 
uvA = xlsread(settings.DisplDataFilePath); 

itMax = settings.itMax;                                                     % maximum number of iterations
tol   = settings.tol;                                                       % tolerance

%% Setup location of inpus files, iteration settings, and initialize variables:                                                                    
mesh = setupmesh;                                                           % setup information

nels = length(mesh.etpl);                                                   % number of elements

error = 2*tol;                                                              % initial error to start the loop
itnum = 0;                                                                  % zero iteration counter

iterationError_E = zeros(itMax,1);
iterationError_v = zeros(itMax,1);
iterationErrorStrainXX = zeros(itMax,1);
iterationErrorStrainYY = zeros(itMax,1);
iterationErrorStrainXY = zeros(itMax,1);
iterationErrorSigXX = zeros(itMax,1);
iterationErrorSigYY = zeros(itMax,1);
iterationErrorSigXY = zeros(itMax,1);
iterationErrorUvXX = zeros(itMax,1);
iterationErrorUvYY = zeros(itMax,1);

%% Iterations:
while error > tol && itnum < itMax                                          % while loop 
    itnum = itnum+1;                                                        % iteration counter
    fprintf('\n%s%8i\n','   iteration number     ',itnum);                  % print iteration number
    [sig,epsH,uv] = LEfe(mesh,itnum);                                       % finite element solver
    [E,v] = updateElasticProp(sig,epsA);                                    % solve for elastic properties 
    mesh.E = E;                                                             % update Young's modulus
    mesh.v = v;                                                             % update Poisson's ratio
    
    error_E = norm(E-YoungRand)/norm(YoungRand);                            % normalised strain error
    iterationError_E(itnum,1) = error_E;
    
    error_v = norm(v-PoissonRand)/norm(PoissonRand);                        % normalised strain error
    iterationError_v(itnum,1) = error_v;
    
    epsH_epsA = epsH-epsA;
    epsH_epsA_XX = epsH_epsA(:,1);
    epsA_XX = epsA(:,1);
    error_strainXX = norm(epsH_epsA_XX)/norm(epsA_XX);
    iterationErrorStrainXX(itnum,1) =  error_strainXX;
    
    epsH_epsA_YY = epsH_epsA(:,2);
    epsA_YY = epsA(:,2);
    error_strainYY = norm(epsH_epsA_YY)/norm(epsA_YY);
    iterationErrorStrainYY(itnum,1) =  error_strainYY;
    
    epsH_epsA_XY = epsH_epsA(:,3);
    epsA_XY = epsA(:,3);
    error_strainXY = norm(epsH_epsA_XY)/norm(epsA_XY);
    iterationErrorStrainXY(itnum,1) =  error_strainXY;                      % save current error to the vector
    
    sig_sigA = sig-sigA;
    sig_sigA_XX = sig_sigA(:,1);
    sigA_XX = sigA(:,1);
    error_sigXX = norm(sig_sigA_XX)/norm(sigA_XX);
    iterationErrorSigXX(itnum,1) =  error_sigXX;
    
    sig_sigA_YY = sig_sigA(:,2);
    sigA_YY = sigA(:,2);
    error_sigYY = norm(sig_sigA_YY)/norm(sigA_YY);
    iterationErrorSigYY(itnum,1) =  error_sigYY;
    
    sig_sigA_XY = sig_sigA(:,3);
    sigA_XY = sigA(:,3);
    error_sigXY = norm(sig_sigA_XY)/norm(sigA_XY);
    iterationErrorSigXY(itnum,1) =  error_sigXY;                            % save current error to the vector
    
    uv_uvA = uv-uvA;
    uv_uvA_XX = uv_uvA(:,1);
    uvA_XX = uvA(:,1);
    error_uvXX = norm(uv_uvA_XX)/norm(uvA_XX);
    iterationErrorUvXX(itnum,1) =  error_uvXX;
    
    uv_uvA_YY = uv_uvA(:,2);
    uvA_YY = uvA(:,2);
    error_uvYY = norm(uv_uvA_YY)/norm(uvA_YY);
    iterationErrorUvYY(itnum,1) =  error_uvYY;
    
end