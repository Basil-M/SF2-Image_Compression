function Y = dwtquant2(Yq, N, dwtstep, rise)

% DWTQUANT2 Reconstruct DWT matrix from quantised values
%
%  The result is the reconstructed values. If rise1 is defined, the first
%  step rises at rise1, otherwise it rises at step/2 to give a uniform
%  quantiser with a step centred on zero.
%  In any case the quantiser is symmetrical about zero.
if size(dwtstep,2) == N + 2
    %rise parameter has been added at the end
    rise = dwtstep(1, N + 2);
    dwtstep = dwtstep(:,1:N+1);
elseif ~exist('rise','var')
    rise = 1;
end

M = size(Yq,1);
Y = zeros(size(Yq));
for i = N:-1:1
    m = M/2^i;
    for k = 1:3
        oy = mod(k, 2);         % will be 1 if k = 1, 3 (so right column)
        ox = (k ~= 1);        % will be 1 if k = 1  (so top row)
        tx = (1+ox*m):(m+ox*m);
        ty = (1+oy*m):(m+oy*m);
        
        Y(tx,ty) = quant2(Yq(tx,ty),dwtstep(k,i),rise*dwtstep(k,i));
    end
end

% Quantise final lowpass image
m = M/2^N;
Y(1:m, 1:m) = quant2(Yq(1:m, 1:m), dwtstep(1, N + 1),dwtstep(1, N + 1)*rise);

