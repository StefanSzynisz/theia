function my_subplot(fig_num,plot_titles,data,nels_x, nels_y)
    
    num_subplots = size(plot_titles,2);                                     % each gap is 20% of each plot
    size_subplot_x = (0.2 + num_subplots + (num_subplots-1)*0.2 + 0.2) * nels_x;
    size_subplot_y = (0.2 + 1.0 + 0.2) * nels_y;
    amplify_fig = min(num_subplots * 1/2,2);                                % amplify the figure for larger number of subplots

    f = figure(fig_num);  close(fig_num); f = figure(fig_num);
    width_fig = f.Position(3); height_fig = f.Position(4);
    scale_fig_size_x = width_fig/size_subplot_x *amplify_fig;
    scale_fig_size_y = height_fig/size_subplot_y *amplify_fig;
%    scale_fig_size = min(scale_fig_size_x,scale_fig_size_y);
%     f.Position(1) = 1/amplify_fig * f.Position(1);
    f.Position(2) = 1/num_subplots * f.Position(2);
    f.Position(3) = size_subplot_x * scale_fig_size_x;
    f.Position(4) = size_subplot_y * scale_fig_size_y;
    
    pos_subplot = zeros(num_subplots,4);                                    % placeholder for the position
    pos_subplot_1 = 0.06;                                                    % positions of the figures
    for ii=1:num_subplots
        pos_subplot(ii,:) = [ pos_subplot_1 0.15 nels_x/size_subplot_x nels_y/size_subplot_y];
        pos_subplot_1 = pos_subplot_1 + (1.2*nels_x)/size_subplot_x;
    end   

    for i=1:num_subplots                                                    % plot each subplot figure
        subplot(1,size(plot_titles,2),i); subplot('Position',pos_subplot(i,:))
        matrix = vector2matrix(data(:,i),nels_y,nels_x);                    % convert vector to matrix for plotting
        imagesc(matrix);colorbar; colormap(hsv); %colormap(hsv);            % plot color map
        axis equal; axis off;
    %     caxis([minSigma maxSigma]);
        title( plot_titles{i} );                                            % add plot title
    %   title( strcat('\',plot_titles{i}) );                                % add plot title
    end
end
