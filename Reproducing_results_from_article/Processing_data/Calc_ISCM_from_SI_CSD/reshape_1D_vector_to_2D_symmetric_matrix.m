
function x = reshape_1D_vector_to_2D_symmetric_matrix(x2)

    % Get size of 2d matrix
    p = 0; k = 0;
    while ~isequal(k,length(x2))
        p = p+1;
        k = (p^2-p)/2;

        if p >1e5
            error('The inputted vector does not correspond to a symmetric 2D matrix, as its length is not equal to (k^2-k)/2, for some integer k (smaller than 100,000)')
        end
    end

    % Build 2D symmetrical matrix
    counter = 1;
    x3 = zeros(p);
    for i = 1:p
        j = p-i;
        x3(:,i) = [zeros(p-j,1) ;x2(counter:counter-1+j)];
        counter = counter + j;
    end
    x4 = x3';
    x = x3 + x4 + eye(p);
end