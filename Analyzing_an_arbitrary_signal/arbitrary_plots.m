  
function [r,f] = arbitrary_plots(t,x,Fs, title_str,thresh_for_corr_inclusion, clip_epsilon)
    % Time plot
    figure,
    hold on
    plot(t,x)
    xlabel('Time (s)')
    ylabel('Amplitude')
    title(['Time-domain ',title_str])
    %plot(t,y)
    
    % Clip
    x = rescale(x,0,1);
    n_noise =   find_matching_end(x,0.01);
    x = x(1:n_noise);

    x = x -mean(x);

%     % Fourier plot
%     figure, line_FFT_xlim(x,Fs,[0 Fs/2]);
%     title(['Frequency-domain ',title_str])

    [S,f,coi] = cwt(x,'morse',Fs);
    figure, cwt(x,'morse',Fs)
    title(['CWT ',title_str])
    freq_limits = [min(f) max(f)];

    
    [r,f,computation_times] = calculate_inter_scale_correlation_matrix(x, Fs, freq_limits, thresh_for_corr_inclusion);
    ylim([min(f) max(f)])
    
    figure,
    plot_2D_colormap(r, f, f, 'log')
    colorbar
    xlabel('Frequency (Hz)')
    ylabel('Frequency (Hz)')
    title(['IFPCM ',title_str])
end