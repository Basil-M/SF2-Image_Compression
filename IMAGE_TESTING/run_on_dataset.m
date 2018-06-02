ref_MSE = @(X) std(X(:) - reshape(quantise(X,17), [numel(X), 1]));
bits = @(X) numel(X)*bpp(X);
%% Laplacian
%Laplacian params
N_lap = 3;
uniform = false; 
h = h_m(2);
c_ratio_py = @(X) py_encbits(X, h, N_lap, findQ(X, ref_MSE(X), N_lap, h, uniform));
%% DCT
%DCT Params
N_dct = 4;
w_subim = 8;
C_dct = dct_ii(N_dct);
Y_dct = @(X) dct_enc(X, N_dct, DCT_q(X,C_dct, ref_MSE(X)));
c_ratio_DCT = @(X) dctbpp(regroup(Y_dct(X),N_dct)/N_dct,w_subim);

%% LBT
%LBT Params
N_lbt = 4;
w_subim = 16;
s = 1.36; %sqrt(2);
C_lbt = dct_ii(N_lbt);

Y_LBT = @(X) dct_enc(X,N_lbt,DCT_q(X,C_lbt,ref_MSE(X),s),s);
c_ratio_LBT = @(X) dctbpp(regroup(Y_LBT(X),N_lbt)/N_lbt, w_subim);
%% DWT
%DWT_params
N_dwt = 7;
Q = @(X) dwt_q(X, ref_MSE(X), N_dwt, false);

c_ratio_DWT = @(X) nlevdwt_bits(X, N_dwt, Q(X));
    
%% All together
c_ratios = @(X) bits(quantise(X,17))./[c_ratio_py(X) c_ratio_DCT(X) c_ratio_LBT(X) c_ratio_DWT(X)];

%% And run

%c_rats = try_on_all(c_ratios);
%save('data_LBT','c_rats')