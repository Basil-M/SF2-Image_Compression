function [Yq, dwtent, nbits] = quantdwt(Y, N, dwtstep)
%QUANTDWT Quantise DWT and calculate entropy/number of bits
%   Y = encoded DWT
%   N = Number of levels
%   dwtstep = Array of quantisation steps

% Quantisation steps top to bottom, left to right so:
%  | ~     k = 1|
%  |k = 2  k = 3|
nu_quant = false;   %Flag - if false, quantise uniformly (no rise)
if length(dwtstep) == N + 2
    %rise parameter has been added at the end
    rise = dwtstep(1, N + 2);
    nu_quant = true;
    dwtstep = dwtstep(:,1:N+1);
end
dwtent = zeros(size(dwtstep));
nbits = 0;

M = size(Y,1);
Yz = zeros(size(Y));
for i = N:-1:1
    m = M/2^i;
    for k = 1:3
        oy = mod(k, 2);         % will be 1 if k = 1, 3 (so right column)
        ox = (k ~= 1);        % will be 1 if k = 1  (so top row)
        tx = (1+ox*m):(m+ox*m);
        ty = (1+oy*m):(m+oy*m);
        
        if nu_quant
            Yq(tx,ty) = quantise(Y(tx,ty),dwtstep(k, i), rise*dwtstep(k,i));
        else
            Yq(tx,ty) = quantise(Y(tx,ty),dwtstep(k, i));
        end
        Yz(tx,ty) = 1;
        dwtent(k,i) = bpp(Yq(tx,ty));
        nbits = nbits + dwtent(k,i)*m^2;
    end
end

% Quantise final lowpass image
m = M/2^N;
Yq(1:m, 1:m) = quantise(Y(1:m, 1:m), dwtstep(1, N + 1));
dwtent(1, N + 1) = bpp(Yq(1:m, 1:m));
nbits = nbits + m^2*dwtent(1, N + 1);
end

