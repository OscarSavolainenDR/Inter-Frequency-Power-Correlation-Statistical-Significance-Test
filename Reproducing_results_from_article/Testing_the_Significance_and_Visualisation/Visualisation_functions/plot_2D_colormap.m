function plot_2D_colormap(matrix, t, f, linear)

    % Function that plots a 2D colormap of a 2D array. Useful for plotting
    % Spectrograms, Wavelet Power Spectrums, and inter-frequency correlation
    % matrices. 
    %
    % Inputs:
    % - matrix: 2D matrix to be plotted.
    % - t: a vector that contains the axis values of the first dimension of the 
    %      2D array. E.g., time stamps, but this can be any arbitrary axis.
    % - f: a vector that contains the axis values of the second dimension of
    %      the 2D array. E.g., frequency axis values, but this can be any
    %      arbitrary axis.
    % - linear: a string, equal to either 'linear', in which case the axes
    %      of the plot are linear, or 'log', in which case both axes are
    %      logarithmic.
    %
    % If "figure" is not called prior to this function, this will
    % overwrite what has been plotted before.
    
    
    %% Tests
    if min(size(t))~=1
        error("'t' should be a 1-D vector")
    end
    
    if min(size(f))~=1
        error("'t' should be a 1-D vector")
    end
    
    if ~isequal(size(matrix),[length(t) length(f)])
        error("'matrix' should be of size '[length(t) length(f)]'")
    end
        
    if ~(strcmp(linear,'linear') || strcmp(linear,'log'))
        error("'linear' variable should be a string equal to 'linear' or 'log'")
    end
    
    %% Plot 2D colormap
    hold on
    surf(f,t,matrix, 'EdgeColor', 'none');
    axis xy; 
    axis tight; 
    colormap(jet); view(0,90);
%     colorbar;
    ax=gca;
    
    %% Set axes
    if strcmp(linear,'linear')
        set(ax,'yscale','linear')
        set(ax,'xscale','linear')
    elseif strcmp(linear,'log')
        set(ax,'yscale','log')
        set(ax,'xscale','log')
    else 
        error("'linear' variable should be a string equal to 'linear' or 'log'")
    end

end
