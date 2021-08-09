function [neural_time, neural_data, behav_time, finger_pos, cursor_pos, target_pos] = Syncing_data(neural_time, neural_data, behav_time, finger_pos, cursor_pos, target_pos, Fs_neural, Fs_behav)
    
    % Trimming data from the start
    starting_temporal_offset = behav_time(1) - neural_time(1);
    if starting_temporal_offset > 0 % Delete some neural data
        offset_samples = round(starting_temporal_offset*Fs_neural);
        neural_time(1:offset_samples) = [];
        neural_data(1:offset_samples,:) = [];
    else % Delete some behavioral data
        offset_samples = round(-starting_temporal_offset*Fs_behav);
        behav_time(1:offset_samples) = [];
        finger_pos(1:offset_samples,:) = [];
        cursor_pos(1:offset_samples,:) = [];
        target_pos(1:offset_samples,:) = [];
    end

    % Trimming data from the end
    end_temporal_offset = behav_time(end) - neural_time(end);
    if end_temporal_offset < 0 % Delete some neural data
        offset_samples = round(-end_temporal_offset*Fs_neural);
        neural_time(end-offset_samples:end) = [];
        neural_data(end-offset_samples:end,:) = [];
    else % Delete some behavioral data
        offset_samples = round(end_temporal_offset*Fs_behav);
        behav_time(end-offset_samples:end) = [];
        finger_pos(end-offset_samples:end,:) = [];
        cursor_pos(end-offset_samples:end,:) = [];
        target_pos(end-offset_samples:end,:) = [];
    end
    
end