%% Uncompress stored neural ISCMs
% There isn't any point saving half of the matrices, as they're
% symmetrical, and the total dataset (mostly null distibutions) takes up 
% 10's of GB. This function takes the compressed matrices and re-formats
% them as individual ISCMs.

function [neural_results] = decompress_neural_ISCMs(stored_neural_ISCMs,neural_parameters,channels_vector)

    neural_results = neural_parameters;

    % Iterate through channels
    clear signal_length
    for channel = channels_vector
        ISCM = reshape_1D_vector_to_2D_symmetric_matrix(stored_neural_ISCMs(:,channel)); 
        neural_results{channel,1}.inter_scale_correlation_matrix = ISCM;
    end
end



