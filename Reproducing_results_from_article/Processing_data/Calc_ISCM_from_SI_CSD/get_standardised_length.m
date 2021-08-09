function [shortened_length] = get_standardised_length(MC_length_vector, minimum_time_gap, session_length)

    % Function that returns the shorter of the 2 closest lengths in the provided vector of
    % possible standardised lengths, for the specified recording.
    % 
    % Inputs:
    % - MC_length_vector: length of standardised lengths, i.e. of the MC
    %   white noise processes.
    % - minimum_time_gap: minimum value that the values in MC_length_vector
    %   need to be seperated by, specified in
    %   "Inter_Scale_Stat_Sig_Sabes_Parameters".
    % - session_length: the length in s of the recording in question
    %
    % Output:
    % - shortened_length: the length in s of the recording in question,
    %   shortened to the closest value in MC_length_vector (given 
    %   minimum_time_gap buffer).
    
    MC_length_vector = sort([MC_length_vector inf]);
    
    for MC_index = 1:numel(MC_length_vector) - 1
        % If between 2 MC_length_vector values, shorten to the smallest one
        if session_length > MC_length_vector(MC_index) -minimum_time_gap && session_length < MC_length_vector(MC_index+1) - minimum_time_gap 
            MC_length(MC_index,1) = 1;
        else 
            MC_length(MC_index,1) = 0;
        end
    end

    % Check the length is long enough to begin with
    if sum(MC_length==1) < 1
        error(['Inputted "session_length" has a length smaller than ',num2str(MC_length_vector(1)-minimum_time_gap),' s.'])
    else
        shortened_length = MC_length_vector(MC_length==1);
    end
end
