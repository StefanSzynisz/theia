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
