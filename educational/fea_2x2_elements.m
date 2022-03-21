% PURPOSE
%     Symbolic solution for 2x2 elements under load, with known BC
%
% DEPENDENCIES:
% 
% RELATED SCRIPTS:
%     stiffness_matrices.m
% Date:
%     03-March-2022
%  ----------------------------------------------------------------
clearvars; close all; clc;
[thisPath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(thisPath); %addpath('functions') 

%%