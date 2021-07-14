function vectSmooth = smoothVector(vect, n_rows, n_cols, filterType, filter_size)
%SMOOTHVECT applies Guassian filter to smooth vector data
%
% Usage: smoothVector(vect,n_rows,n_cols, filter_size)
% vect - input vector
% n_rows - number of rows in the grid of elements
% n_cols - number of columns in the grid of elements
% filterType = 'Gaussian' or 'MovingAverage'
% filter_size = 3, 5, or any higher odd number. 

    
    matrix = vector2matrix(vect,n_rows,n_cols); % first vector to process

    % Moving average:
    % https://uk.mathworks.com/matlabcentral/answers/358411-i-need-a-code-that-produce-a-moving-average-matrix-with-a-5-5-window
    % K = 1/9*ones(3); % K = 1/16*ones(4); K = 1/25*ones(5); % moving average
    % Guassian blur 3 x 3:
    %--------------

    if strcmp(filterType,'MovingAverage')
        K = ones(filter_size) ./ (filter_size)^2;
        % For example:
        %K = 1/9 * ...
        %    [1 1 1;
        %    1 1 1;
        %    1 1 1]; 
        matrixSmoothTMP = conv2(matrix,K,'same');   % Apply the filter

    elseif strcmp(filterType,'Gaussian')
        sigma = 1.0; %such that the filter is spread across the entire matrix
        % Filter can be also generated using fspecial and applied using conv2
        %K = fspecial('gaussian', filter_size, sigma);
        %matrixSmoothTMP = conv2(matrix,K,'same');   % Apply the filter
    
        % Gaussian blur 3 x 3 (an example of the weights):
        %K = 1/16 * ...
        %    [1 2 1;
        %    2 4 2;
        %    1 2 1];
        % Apply the filter using built-in function:
        matrixSmoothTMP = imgaussfilt(matrix,sigma,'FilterSize',filter_size);
    end

    % Remove smoothed results from edge_width pixels around the frame
    % These edge results do not seem right after the smoothing.
    % Check if the filter size is odd:
    if rem(filter_size,2)  % reminder after division by 2 is 1 = true
        edge_width = (filter_size - 1)/2;
        % fprintf('%i cell(s) around the edge are beyond the reach of the smoothing algorithm.\n', edge_width);
    else
        fprintf('Only odd numbers are allowed as filter_size and %i is not an odd number.\n\n', filter_size);
    end

    %Keep the original edges of the images:
    matrixSmooth = keepOriginalEdges(matrixSmoothTMP,matrix,edge_width);
    % matrixSmooth = matrixSmoothTMP;

    vectSmooth = matrix2vector(matrixSmooth,n_rows,n_cols);


    % FUNCTIONS 
    % =========
    function matrix = vector2matrix(vector,n_rows,n_cols)  
        matrix = zeros(n_rows, n_cols); % Pre-allocate the matrix
        for k = 1:1:n_rows
            matr_row = n_rows+1-k;
            vect_indx = 1 + (k-1)*n_cols;
            vect_indx_end = vect_indx + (n_cols -1);
            matrix(matr_row,:) = vector(vect_indx:vect_indx_end);
        end
    end

    function matrixOUT = keepOriginalEdges(matrixSmooth, matrixOrig,edge_width)
        matrixOUT = matrixSmooth;
        % but keep the first 'edge_width' rows, last two rows, first two columns and
        % last two columns
        
        edge_subt = edge_width-1;
        matrixOUT(1:edge_width,:) = matrixOrig(1:edge_width,:);
        matrixOUT(end-edge_subt:end,:) = matrixOrig(end-edge_subt:end,:);
        matrixOUT(:,1:edge_width) = matrixOrig(:,1:edge_width);
        matrixOUT(:,end-edge_subt:end) = matrixOrig(:,end-edge_subt:end);
    end

    function vector = matrix2vector(matrix,n_rows,n_cols)  
        vector = zeros(n_rows*n_cols,1);
        for k = 1:1:n_rows
            matr_row = n_rows+1-k;
            vect_indx = 1 + (k-1)*n_cols;
            vect_indx_end = vect_indx + (n_cols -1);
    
            vector(vect_indx:vect_indx_end) = matrix(matr_row,:);
        end
    end

end