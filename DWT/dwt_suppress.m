function [Y_sup] = dwt_suppress(Y, N, N_sup)
%DWT_SUPPRESS Suppress N_supress sub-images DWT output Y 
%   Detailed explanation goes here
rats = dwt_q_ratios(size(Y), N);        %Suppress in increasing order of energy
rats(2:3, N+1) = -1;
% functions to get range of pixels corresponding to index in rats
k = 3; i = N;
M = size(Y,1);

for n = 1:N_sup
    m = M/2^i;
    oy = mod(k, 2);         % will be 1 if k = 1, 3 (so right column)
    ox = (k ~= 1);        % will be 1 if k = 1  (so top row)
    tx = (1+ox*m):(m+ox*m);
    ty = (1+oy*m):(m+oy*m);
    
    Y(tx,ty) = quantise(Y(tx,ty), 10, 10000);
    rats(k, i) = -1;
    
    if k == 1
        k = 3;
        i = i - 1;
    else
        k = k - 1;
    end 
end
Y_sup = Y;
end

