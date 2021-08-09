function [x2] = reshape_2D_symmetric_matrix_to_1D(x)
    x2 = [];
    for i = 1:length(x)
        x2 = [x2; x(i+1:end,i)];
    end
end
