function plot_2D_colormap_nan_ignore(matrix, t, f, linear, colors, colorbar_bool, plot_abs_corr, caxes)

    % Function that plots a 2D colormap of a 2D array. Useful for plotting
    % Spectrograms, Wavelet Power Spectrums, and inter-frequency correlation
    % matrices. Relative to plot_Spectrogram, this has the advantage of NaN 
    % values not graphically over-writing non-NaN values. Graphically, this
    % means that all significant values get shown.
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
    % - colors: MATLAB standard colormaps, e.g. jet, cool, etc.
    % - colorbar_bool: a boolean that defines whether we plot the colorbar or not. Sometimes we don't want to,
    %      e.g. if we want to plot a universal colorbar across different subplots.
    % - plot_abs_corr:
    % - caxes:
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
    
    %% Creating custom colormap, with NaN values = black
    gg = colors;
    clear cc_temp_triple cc_temp jjj
    matrix(1,1) = caxes(1); matrix(2,2) = caxes(2);  % add global min and max to the matrix, so the re-scaling is correct
    cc_temp = round(rescale(matrix,1,256));
    cc_temp(1,1) = NaN; cc_temp(2,2) = NaN;
    for hh = 1:length(t)
        for hhh = 1:length(f)
            jjj = cc_temp(hh,hhh);
            if isnan(jjj) % if NaN, set color as equal to black
                cc_temp_triple(hh,hhh,:) = [0 0 0];
            else          % otherwise, usual color
                cc_temp_triple(hh,hhh,:) = gg(jjj,:);
            end
        end
    end
    
    %% Plot 2D colormap
    surf(f,t,ones(length(f),length(f)),cc_temp_triple, 'EdgeColor', 'none');
    axis xy; 
    axis tight; 
    colormap(jet); view(0,90);
    if colorbar_bool % true = we plot the colorbar
        hc1 = colorbar;
        if ~plot_abs_corr
            set(get(hc1,'label'),'string','Correlation');
        else
            set(get(hc1,'label'),'string','| Correlation |');
        end
        caxis([min(min(matrix)) max(max(matrix))])
    end

    %% Set axes
    ax = gca;
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
