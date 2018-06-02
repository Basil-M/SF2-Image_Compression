function [dwtstep] = dwt_q_ratios_X(X, N)
%DWT_Q_RATIOS Finds equal MSE Q ratios for n level dwt
%   sz = Size of image
%   N = number of levels of dwt ratio

%Make fake Y
Y = nlevdwt(X, N);
sz = size(Y);

E = @(Z) sum(Z(:).^2);  %Calculate energy
E_ki = Inf*ones(3, N + 1);

%Coords of impulse are the same each time just scaled by 2
co = zeros(3,4);
co(1, :) = [1,0.5*sz(1),0.5*sz(2)+1, sz(2)];    % k = 1, top right corner
co(2, :) = [0.5*sz(1)+1,sz(1), 1,0.5*sz(2)];    % k = 2, bottom left corner
co(3, :) = [0.5*sz(1)+1,sz(1), 0.5*sz(2)+1,sz(2)];    % k = 3, bottom right corner

for i = 1:N
    for k = 1:3
        E_ki(k, i) = E(Y(co(k,1):co(k,2), co(k,3):co(k,4)))/numel(Y(co(k,1):co(k,2), co(k,3):co(k,4)));
    end
    co = co/2;
end

%final lowpass
m = sz(1)/2^N;
E_ki(1, N + 1) = E(Y(1:m,1:m))/m^2;
dwtstep = 1./sqrt(E_ki);
dwtstep = dwtstep/max(dwtstep(:)); %normalise
end

