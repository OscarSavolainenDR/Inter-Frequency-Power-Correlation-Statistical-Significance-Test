function [map_sig] = return_FDR_stat_map(h,map2,f,stat,plot_figure)

    % Function that returns the 2D map from 1D vector
    % 
    % Inputs:
    % - h: vector of rejected hypotheses
    % - map2: 2D matrix, where the value indicates what hypothesis it
    %   corresponds to.
    % - f: vector of frequencies in the inter-frequency correlation matrix
    % - stat: there are 10 tests in the normlaity testing, stat is the
    %   index of the current tests, e.g. 1-10.
    % - plot_figure: boolean value, whether to plot the map of
    %   not-statistically non-normal null distributions.
    %
    % Outputs: 
    % - mag_sig: map of not-statistically non-normal null distributions.
            
    temp = find(h==1);
    temp = temp(temp>0);
    
    %% Get 2D map of hypotheses
    map_sig = zeros(length(f));
    counter = 0;
    for k1 = 1:length(f)
        for kk = 1:k1-1
            counter = counter + 1;
            map(k1,kk) = counter;
        end
    end

    %% Find which hypotheses in the map are statistically significant
    for i = 1:length(temp)
        [row, col] = find(map==map2(temp(i)));
        Row(i,1) = row; Col(i,1) = col;
    end

    %% Create new map only of statisically significant hypotheses
    for i = 1:length(Row)
        map_sig(Row(i),Col(i)) = 1;
        map_sig(Col(i),Row(i)) = 1;
    end
    
    %% Plot statistically significant inter-scale correlation map
    if plot_figure
        plot_2D_colormap(map_sig,f,f,'log')
        xlabel('Frequency (Hz)'); ylabel('Frequency (Hz)')
        title(['FDR indicator map, stat test nb. ',num2str(stat)])
    end

end
    
