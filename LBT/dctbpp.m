function [n_bits] = dctbpp(Y,N)
%dct_ENERGIES Calculates bits per pixel (entropy) of dct with N x N filter
%   Y = transformed image
%   N = size of transform
n_bits = 0;
E = @(X) sum(X(:).^2);

l_r = size(Y,1)/N; l_c = size(Y,2)/N;

for r = 1:N
    r_lst = ((r - 1)*l_r + 1):r*l_r;        % Row indices
    for c = 1:N
        c_lst = ((c - 1)*l_c + 1):c*l_c;    % Column indices
        n_bits = n_bits + bpp(Y(r_lst, c_lst))*l_c*l_r;
    end
end

