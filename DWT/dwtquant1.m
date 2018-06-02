function Yq = dwtquant1(Y, N, dwtstep, rise)

%DWTQUANT1 Quantise DWT and calculate entropy/number of bits
%   Y = encoded DWT
%   N = Number of levels
%   dwtstep = Array of quantisation steps

%  The result is the quantised integers Q. If rise1 is defined,
%  the first step rises at rise1, otherwise it rises at step/2 to
%  give a uniform quantiser with a step centred on zero.
%  In any case the quantiser is symmetrical about zero.


% Quantisation steps top to bottom, left to right so:
%  | ~     k = 1|
%  |k = 2  k = 3|
if length(dwtstep) == N + 2
    %rise parameter has been added at the end
    rise = dwtstep(1, N + 2);
    dwtstep = dwtstep(:,1:N+1);
elseif ~exist('rise','var')
    rise = 1;
end

M = size(Y,1);
Yz = zeros(size(Y));
for i = N:-1:1
    m = M/2^i;
    for k = 1:3
        oy = mod(k, 2);         % will be 1 if k = 1, 3 (so right column)
        ox = (k ~= 1);        % will be 1 if k = 1  (so top row)
        tx = (1+ox*m):(m+ox*m);
        ty = (1+oy*m):(m+oy*m);
        
        Yq(tx,ty) = quant1(Y(tx,ty),dwtstep(k,i),rise*dwtstep(k,i));
    end
end

% Quantise final lowpass image
m = M/2^N;
Yq(1:m, 1:m) = quant1(Y(1:m, 1:m), dwtstep(1, N + 1),dwtstep(1, N + 1)*rise);

