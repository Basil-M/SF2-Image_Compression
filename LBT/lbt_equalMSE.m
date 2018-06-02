% N is the number of layers, i is the size of the DCT transform
function Q = lbt_equalMSE(W, N_LEVELS, N, s)
Q = zeros([N_LEVELS, N, N]);
Y_zero = zeros(W);
E = @(Z) sum(Z(:).^2);  %Calculate energy


for n = 1:N_LEVELS
    for r = 1:N
        for c = 1:N
            Y = Y_zero;
            x_rc = r + N*ceil(0.5*(W-r)/N);
            y_rc = c + N*ceil(0.5*(W-c)/N);
            Y(x_rc,y_rc) = 100;
            X_rec = lbt_dec(Y, N, s);
            Q(n,r,c) = 100^2/E(X_rec);
        end
    end
end

