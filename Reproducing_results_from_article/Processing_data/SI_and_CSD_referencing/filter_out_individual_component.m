function Y = filter_out_individual_component(Y, Fs, component_freq, component_neigh)

    % Spectral interpolation and phase randomisation of single component.
    % Inputs:
    % - Y: Fourier transform of data, with components to be interpolated across
    % - Fs_neural: sampling frequency of the data (generally neural data)
    % - component_freq: cell matrix that contains lower and upper limits of the components to be interpolated across.
    % - component_neigh: scalar value, how many frequencies to use in calculating the average value on either side of the interpolation. Practically, this value is unimportant, as the data is later pre-whitened so the amplitude of the interpolation is insignificant.
    % 
    % Output:
    % - Y: Fourier transform of data, with specified components removed


    % Tests:
    if component_neigh < 0
        error('component_neigh must be a positive value. \n');
    end
    if Fs < 0
        error('Fs must be a positive value. \n');
    end
    if rem(numel(component_freq),2) ~= 0
        error('component_freq should have 2k elements, where k is the number of components to interpolate .\n');
    end
    if length(Y(:,1)) < length(Y(1,:)) % if the matrix is entered the wrong way, i.e. the number of frequencies is smaller than the number of channels
        Y = Y';
    end

    % Get frequency vector
    L = length(Y(:,1));
    nb_channels = length(Y(1,:));
    frq = Fs*linspace(0,1,L);

    % Iterate through components, to interpolate over them
    counter = 0;
    component_nb = length(component_freq)/2; % nb of components to interpolate across
    while counter < component_nb
        counter = counter + 1;
        
        % Frequencies to interpolate
        f2int = [component_freq((counter-1)*2+1) component_freq(counter*2)];
    
        % Frequencies to use in the interpolation average
        f4int = [f2int(1)-component_neigh f2int f2int(2)+component_neigh];

        % Find indices of frequencies to interpolate
        [~,gg] = min(abs(frq-f2int(1))); % find closest measured frequency to component (left of component)
        [~,gg2] = min(abs(frq-f2int(2))); % right of component)
        smpl2int = [gg gg2]; % range of frequencies that will be interpolated

        % Find indices of frequencies used for interpolation
        [~,gg3] = min(abs(frq-f4int(1))); % find closest measured frequency to theoretical line-noise harmonic
        [~,gg4] = min(abs(frq-f4int(4)));
        smpl4int = [gg3 smpl2int gg4]; % range of frequencies that will be used for interpolation

        for channel = 1:nb_channels
            % To add random phase shifts (negative for conjugates), preserve DC
            % offset. Here we calculate random phases.
            rnd_theta = -pi + (2*pi).*rand(1,smpl2int(2)-smpl2int(1)+1);

            P2 = abs(Y(:,channel)); % FFT abs

            % Interpolation magnitude
            x = [mean(P2(smpl4int(1):smpl4int(2))) mean(P2(smpl4int(3):smpl4int(4)))]; % 2 meaned values of either side of interpolation
            coefficients = polyfit([smpl2int(1), smpl2int(2)], [x(1), x(2)], 1); % 1st order polynomial to span the interpolation "gap"
            temp = coefficients(1)*[smpl2int(1):smpl2int(2)] + coefficients(2); % a*x + b form, where x is the harmonic gap to interpolate over

            % Patch spectra with interpolation
            Y(smpl2int(1):smpl2int(2),channel) = temp.*exp(1i*rnd_theta); % interpolated magnitude with randomised phase
            Y(L-smpl2int(2):L-smpl2int(1),channel) = flip(temp.*exp(-1i*rnd_theta)); % same for negative frequencies
        end
    end

end
