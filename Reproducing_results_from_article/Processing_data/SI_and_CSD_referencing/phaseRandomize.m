function randX = phaseRandomize(X)

% Returns a phase-randomized version of the input data X. If X is an array,
% each row is treated as an independant time series, and columns represent
% sample points. 
[N,L]=size(X);
Y=fft(X,[],2); % Get spectrum

% Add random phase shifts (negative for conjugates), preserve DC offset
rnd_theta= -pi + (2*pi).*rand(N,floor(L/2-1)); 
Y(:,2:L/2)=Y(:,2:L/2).*exp(1i*rnd_theta);
Y(:,L/2+2:L)=Y(:,L/2+2:L).*exp(-1i*flip(rnd_theta,2));
% return phase-randomized data
randX =ifft(Y,[],2);
end