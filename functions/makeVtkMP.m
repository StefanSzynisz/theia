function makeVtkMP(mpC,sig,strain,uvw,E,v,mpFileName)

%VTK output file generation: Gauss point data
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   15/01/2019
% Description:
% Function to generate a VTK file containing the Gauss point data
%
%--------------------------------------------------------------------------
% MAKEVTKMP(mpC,sig,strain,uvw,E,v,mpFileName)
%--------------------------------------------------------------------------
% Input(s):
% mpC        - point coordinates (nmp,nD)
% sig        - stresses (nmp,*)
% strain     - strains (nmp,*)
% uvw        - displacements (nmp,nD)
% E          - Young's modulus (nmp,1)
% v          - Poisson's ratio (nmp,1)
% mpFileName - VTK file name, for example 'mpData.vtk'  
%--------------------------------------------------------------------------
% See also:
% 
%--------------------------------------------------------------------------

[nmp,nD]=size(mpC);                                                         % no. points and dimensions

fid=fopen(mpFileName,'wt');
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'MATLAB generated vtk file, WMC\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %i double\n',nmp);

%% position output 
if nD<3
    mpC = [mpC zeros(nmp,3-nD)];  
end
fprintf(fid,'%f %f %f \n',mpC');
fprintf(fid,'\n');

fprintf(fid,'POINT_DATA %i\n',nmp);
%% stress output
if nD==3
    fprintf(fid,'SCALARS sigma_xx FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',sig(:,1));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS sigma_yy FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',sig(:,2));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS sigma_zz FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',sig(:,3));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS sigma_xy FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',sig(:,4));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS sigma_yz FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',sig(:,5));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS sigma_zx FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',sig(:,6));
    fprintf(fid,'\n');
elseif nD==2
    fprintf(fid,'SCALARS sigma_xx FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',sig(:,1));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS sigma_yy FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',sig(:,2));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS sigma_xy FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',sig(:,3));
    fprintf(fid,'\n');
elseif nD==1
    fprintf(fid,'SCALARS sigma_xx FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',sig(:,1));
    fprintf(fid,'\n');
end

%% strain output
if nD==3
    fprintf(fid,'SCALARS eps_xx FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',strain(:,1));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS eps_yy FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',strain(:,2));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS eps_zz FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',strain(:,3));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS eps_xy FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',strain(:,4));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS eps_yz FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',strain(:,5));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS eps_zx FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',strain(:,6));
    fprintf(fid,'\n');
elseif nD==2
    fprintf(fid,'SCALARS eps_xx FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',strain(:,1));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS eps_yy FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',strain(:,2));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS eps_xy FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',strain(:,3));
    fprintf(fid,'\n');
elseif nD==1
    fprintf(fid,'SCALARS eps_xx FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',strain(:,1));
    fprintf(fid,'\n');
end


%% displacement output
if nD==3
    fprintf(fid,'SCALARS u_x FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw(:,1));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS u_y FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw(:,2));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS u_z FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw(:,3));
    fprintf(fid,'\n');
elseif nD==2
    fprintf(fid,'SCALARS u_x FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw(:,1));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS u_y FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw(:,2));
    fprintf(fid,'\n');
elseif nD==1
    fprintf(fid,'SCALARS u_x FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw);
    fprintf(fid,'\n');
end

%% material properties output
fprintf(fid,'SCALARS E FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',E);
fprintf(fid,'\n');

fprintf(fid,'SCALARS v FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',v);
fprintf(fid,'\n');

fclose('all');