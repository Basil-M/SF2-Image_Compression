function [dwtstep] = dwt_q_ratios(sz, N)
%DWT_Q_RATIOS Finds equal MSE Q ratios for n level dwt
%   sz = Size of image
%   N = number of levels of dwt ratio

%Make fake Y
Y_zeros = zeros(sz);

%Set impulse value
E_imp = 100;

E = @(Z) sum(Z(:).^2);  %Calculate energy
E_ki = zeros(3, N + 1);

%Coords of impulse are the same each time just scaled by 2
coords = zeros(3,2);
coords(1, :) = [0.25*sz(1), 0.75*sz(2)];    % k = 1, top right corner
coords(2, :) = [0.75*sz(1), 0.25*sz(2)];    % k = 2, bottom left corner
coords(3, :) = [0.75*sz(1), 0.75*sz(2)];    % k = 3, bottom right corner

for i = 1:N
    for k = 1:3
        Y = Y_zeros;
        Y(coords(k,1), coords(k,2)) = E_imp;
        X_t = nlevidwt(Y, N);
        E_ki(k, i) = E(X_t)/E_imp^2;
    end
    coords = coords/2;
end

%final lowpass
Y = Y_zeros;
m = sz(1)/2^(N+1);
Y(m,m) = E_imp;
X_t = nlevidwt(Y, N);
E_ki(1, N + 1) = E(X_t)/E_imp^2;
dwtstep = 1./sqrt(E_ki);
dwtstep = dwtstep/dwtstep(1,1); %normalise
%E_ki(1, N + 1) = 0;
end

