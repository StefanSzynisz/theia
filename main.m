%plane stress elastic property solver 
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   15/04/2020
% Description:
% Script to determine the elastic properties of a material based on a known
% strain field, epsA, and predicted stress field from linear elastic finite
% element analysis.  The script iterates to find the elastic properties of
% the material until it converges within a given tolerance or reaches a
% maximum number of iteraitons (default set at 50). 
%
% The script requires a xlsx file with the known strain information:
% Strain.xlsx
%
% The physical problem should be defined in setupmesh.m 
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
% SETUPMESH         - problem set up 
% LEFE              - finite element solver
% UPDATEELASTICPROP - elastic property update
%--------------------------------------------------------------------------
clearvars; close all; clc;
[thisPath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(thisPath);

clear; tic;
addpath('functions');                                                       % add path to functions folder
mesh = setupmesh;                                                           % setup information

nels = length(mesh.etpl);                                                   % number of elements
epsA = xlsread('Strain.xlsx');                                              % read strain data

itMax = 50;                                                                  % maximum number of iterations
tol   = 1e-9;                                                               % tolerance
error = 2*tol;                                                              % initial error
itnum = 0;                                                                  % zero iteration counter

while error > tol && itnum < itMax                                          % while loop 
    itnum = itnum+1;                                                        % iteration counter
    fprintf('\n%s%8i\n','   iteration number     ',itnum);                  % print iteration number
    [sig,epsH] = LEfe(mesh,itnum);                                          % finite element solver
    [E,v] = updateElasticProp(sig,epsA);                                    % solve for elastic properties
    mesh.E = E;                                                             % update Young's modulus
    mesh.v = v;                                                             % update Poisson's ratio
    error = norm(epsH-epsA)/norm(epsA);                                     % normalised strain error
    fprintf('%s%8.3e\n','   error                   ',error);               % return strain error
end