function [P2,f]= line_FFT_xlim(signal,Fs,xl)

    % Plots FFT
    % Inputs:
    % - signal: time-series for which the FFT is to be plotted.
    % - Fs: sampling frequency of signal.
    % - xl: 2-element vector, acts as xlim([xl(1) xl(2)]) for the Fourier 
    %   domain plot. This limits the size of the plots to speed up the
    %   computation.
    %
    % Outputs:
    % - P2: returns the Fourier spectra.
    % - f:  the frequency axis vector.
    %
    % If figure is not called prior to this function, it will over-write
    % what has been plotted before.
    
    %% Calculate FFT
    L = length(signal); % Length of signal
    Y = fft(signal); % FFT
    P2 = abs(Y/L); 
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    
    %% Plot FFT
    [~,lower_limit] = min( abs( f-min(xl) ) );
    [~,upper_limit] = min( abs( f-max(xl) ) );
    line(f(lower_limit:upper_limit),P1(lower_limit:upper_limit)) 
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')
    
end
