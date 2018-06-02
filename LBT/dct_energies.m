function [en] = dct_energies(Y,N)
%dct_ENERGIES Calculates energies of dct 
%   Y = encoded dct
%   N = size of filter
en = zeros([N,N]);
E = @(X) sum(X(:).^2);
l_r = size(Y,1)/N; l_c = size(Y,2)/N;
Y = regroup(Y,N)/N;
for r = 1:N
    r_lst = ((r - 1)*l_r + 1):r*l_r;        % Row indices
    for c = 1:N
        c_lst = ((c - 1)*l_c + 1):c*l_c;    % Column indices
        en(r,c) = E(Y(r_lst, c_lst));
    end
end

end


