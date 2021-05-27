% This is a function which can directly read data from processed ncorr data
% Ncorr is a open source MATLAB platform that does digital image
% correlation from images of any source. More info can be found here: http://www.ncorr.com/
% In this function, a 

%      Y
%      ^
%      |                                |
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
%                            (nx * ny)




function [epsA,nx,ny,subsetspace] = ncorr_read(filename)
% input filename, output ordered strain components
% convention, e11,e22,e12
% output is a lxnx3 matrix where l is the number of frames n is number of
% pixels 
% nx, ny is the number of elements in x and y direction
% if dimension is specified in ncorr, subsetspace is the size of the
% elemetns in which ever dimension is specified
    load(filename);

    strain_data=data_dic_save.strains;
    strain_roi = strain_data.plot_exx_ref_formatted;
    subsetspace=(data_dic_save.dispinfo.spacing+1)*data_dic_save.dispinfo.pixtounits;
    
    % extract the rectangular bound of the space
    [rid,cid] = find(strain_roi);
    bound = [min(rid),max(rid),min(cid),max(cid)];
    nx = bound(4)-bound(3)+1;
    ny = bound(2)-bound(1)+1;
    
    l=length(current_save);
    
    epsA = zeros(l,nx*ny,3);
    
    for i =1:l
        % read frame strain data
        strain = strain_data(i);
        exx = strain.plot_exx_ref_formatted;
        eyy = strain.plot_eyy_ref_formatted;
        exy = strain.plot_exy_ref_formatted;
        exx=exx(bound(1):bound(2),bound(3):bound(4));
        eyy=eyy(bound(1):bound(2),bound(3):bound(4));
        exy=exy(bound(1):bound(2),bound(3):bound(4));
        
        % reformat strain according to convention
        exx = reshape(transpose(flip(exx,2)),[],1);
        eyy = reshape(transpose(flip(eyy,2)),[],1);
        exy = reshape(transpose(flip(exy,2)),[],1);
        
        epsA(i,:,:) = [exx,eyy,exy];
        
    end
end
