%% experimenting with different schemes for final scheme
%% Optimal parameters
% DCT
N_dct = 8; N_enc_dct = 8;
% LBT
N_lbt = 4; s_lbt = sqrt(2); N_enc_lbt = 16;
%DWT
N_dwt = 6; constStep_dwt = false;

%% multi level LBT 
N=4; 
s = sqrt(2);
%ENCODE
Y = lbt_justenc(X, N, s);


W = 256;
H = 256;

kmax = W/N - 1; % the number of components to pick out in the width

lmax = H/N - 1; % the number of components to pick out in the height

A = 0;

% pick out the DC component
A = Y(1:4:256, 1:4:256)/N;

% encode the dc image
B = lbt_justenc(A, N*2, s);


Y(1:4:256, 1:4:256) = B;




Y(1:4:256, 1:4:256) = lbt_dec(Y(1:4:256, 1:4:256), N*2, s)*N;

Z = lbt_dec(Y, N, s);

draw(Z);
std(X(:)-Z(:))


%%
size(A)
figure(1)
draw(A);

Yr = regroup(Y, N)/N;
figure(2);
draw(Yr);
Yrdc = Yr(1:lmax+1, 1:kmax+1);
draw(Yrdc);
rmsError = std(Yrdc(:)-A(:))

%% 3 level lbt
% encode
load lighthouse;
N = 4;
s = sqrt(2);

Y = lbt_justenc(X-128, N, s);
% performing a 2 level lbt

A = 0;

% pick out the DC component
A = Y(1:4:256, 1:4:256)/N;

    % encode the dc coefficients (with a 2Nx2N dct block if you put 2*N
    % into lbt_dec)
    B = lbt_justenc(A, N, s);

    % put in if statement here so can make the third level optional.
    % pick out dc coefficients of 2nd level
    C = B(1:4:64, 1:4:64)/N;

    % encode the dc coefficients 
    D = lbt_justenc(C, N, s);

    B(1:4:64,1:4:64) = C;


Y(1:4:256, 1:4:256) = B;

% decode 
Zi = Y;
% add in if statement for this
    % decode the encoded 2nd level dc coefficients
    Zj = Zi(1:4:256, 1:4:256);
    Zj(1:4:64,1:4:64) = lbt_dec(Zj(1:4:64,1:4:64),N,s)*N;
    Zi(1:4:256, 1:4:256) = Zj;
% decode the encoded dc coefficients (with a 2Nx2N dct block if you put 2*N
% into lbt_dec)
Zi(1:4:256, 1:4:256) = lbt_dec(Zi(1:4:256, 1:4:256), N, s)*N;

Z = lbt_dec(Zi, N, s);

