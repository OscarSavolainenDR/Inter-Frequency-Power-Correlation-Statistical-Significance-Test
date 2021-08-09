%% Compress 2D neural ISCMs
% There isn't any point saving half of the matrices, as they're
% symmetrical, and the total dataset (mostly null distibutions) takes up 
% 10's of GB.

function [stored_neural_ISCMs] = compress_neural_ISCMs(channels_vector,neural_results,p)

%     % Iterate through sessions
%     for session = sessions_vector

        % Initialise 2D matrix, where each column represents one channel's ISCM in
        % the same session
        stored_neural_ISCMs = zeros((p^2-p)/2,length(channels_vector));

        % Iterate through channels
        clear signal_length stored_neural_ISCMs
        for channel = channels_vector
            ISCM = neural_results{channel, 1}.inter_scale_correlation_matrix;
            [x2] = reshape_2D_symmetric_matrix_to_1D(ISCM);
            if isequal(length(x2),(p^2-p)/2)
                stored_neural_ISCMs(:,channel) = x2;
            else
                stored_neural_ISCMs(:,channel) = zeros((p^2-p)/2,1);
            end
        end
%     end
end



