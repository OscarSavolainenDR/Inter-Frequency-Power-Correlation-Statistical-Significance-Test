function index = find_matching_end(x,epsilon)

    % Finds matching end of sequence, so that the beginning is within some
    % value epsilon of the end. Returns the index value of the end that
    % matches the beginning. Used for clipping prior to pre-whitening to
    % prevent edge effects due to cyclic discontinuities.

    x_beginning = abs(x - x(1));
    temp = single(x_beginning<epsilon);
    [~,index]= max(cumsum(temp));
    
end
