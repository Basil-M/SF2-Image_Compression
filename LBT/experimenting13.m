
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
%draw(Yrdc);
rmsError = std(Yrdc(:)-A(:))

%% 3 level lbt
% encode
load lighthouse;
N = 4;
s = sqrt(2);

Y = lbt_justenc(X-128, N, s);
% performing a 3 level lbt

A = 0;

% pick out the DC component
A = Y(1:4:256, 1:4:256)/N;

    % encode the dc coefficients (with a 2Nx2N dct block if you put 2*N
    % into lbt_dec)
    B = lbt_justenc(A, N, s);

    % put in if statement here so can make the third level optional.
    % pick out dc coefficients of 2nd level
    C = B(1:4:64, 1:4:64)/N;
    draw(C);
    % encode the dc coefficients 
    D = lbt_justenc(C, N, s);

    
    %{
    C1 = lbt_dec(D,N,s);
    B1(1:4:64, 1:4:64) = C1*N;
    A1 = lbt_dec(B, N, s);
    Y(1:4:256,1:4:256) = A1*N;
    X1 = lbt_dec(Y, N, s);
    draw(X1);
    std(X(:)-X1(:))
    %}

    B(1:4:64,1:4:64) = D;


Y(1:4:256, 1:4:256) = B;

% decode 
Zi = Y;
% add in if statement for this
    % decode the encoded 2nd level dc coefficients
    Zj = Zi(1:4:256, 1:4:256);
    
    %Zj(1:4:64,1:4:64) = lbt_dec(Zj(1:4:64,1:4:64),N,s)*N;
    Zi(1:16:256,1:16:256)=lbt_dec(Zi(1:16:256,1:16:256),N,s)*N;
    %Zi(1:4:256, 1:4:256) = Zj;
% decode the encoded dc coefficients (with a 2Nx2N dct block if you put 2*N
% into lbt_dec)
Zi(1:4:256, 1:4:256) = lbt_dec(Zi(1:4:256, 1:4:256), N, s)*N;

Z = lbt_dec(Zi, N, s);

%% different levels of ratio between the quantisation of the first level of the lbt and the 2nd level of the lbt.
%% currently no quantisation is being used for the 3rd level of the lbt
N = 4;
M= 16;
rise1 = 1;
s = sqrt(2);


i=1;
ssimval1=0;
ssimval2=0;
for ratio=0.05:0.05:0.5
    load lighthouse;
    [ssimval1(i), ~, Z1, ~] = jpegencdeclbtnlev(X, N, M, rise1, s, true, 16, ratio);
    load bridge;
    [ssimval2(i), ~, Z2, ~] = jpegencdeclbtnlev(X, N, M, rise1, s, true, 16, ratio);
    i = i+1;
end



%% This didn't work 
evalc(' ssimval = @(ratio)jpegencdeclbtnlev(X, N, M, rise1, s, true, 16, ratio); ratio_opt = fminsearch(ssimval, 0.2);');
%%
plot(0.05:0.05:0.5, ssimval);
xlabel('Ratio of quantisation'); % for the 1st level of lbt to the 2nd level of lbt
ylabel('ssim value');

%% Using the bridge image and taking the ratio which gave the optimum value (0.3). Testing the best ratio for the 3rd level

load bridge;
N = 4;
M= 16;
rise1 = 1;
s = sqrt(2);
ratio = 0.3;
i=1;
ssimval =0;
for ratio2=0:0.05:0.5
    [ssimval2(i), rmsError2, Z2, q_opt2] = jpegencdeclbtnlev(X, N, M, rise1, s, true, 16, ratio, ratio2);
    i = i+1;
end

%% 
plot(0:0.05:0.5, ssimval);
xlabel('Ratio of quantisation'); % for the 1st level of lbt to the 2nd level of lbt
ylabel('ssim value');

%% Equal energy for one level of lbt - basing the quantisation on the image itself
% for a lbt using a 4x4 dct transform, there will be 16 subimages.
N=4;

energy = 0;

for i=1:4
    for j=1:4
        Xsub = X(i:N:256,j:N:256);
        energy(i,j) = sum(Xsub(:).^2);
    end
end 
energy = energy/sum(energy(:));

en = dct_energies(X, N);
en = en/sum(en(:)); 

std(en(:)-energy(:))
%}
%% 
N = 4;
s = sqrt(2);
step = 17;
rise1 = 1;






Y = lbt_justenc(X, N, s);

en = dct_energies(Y,N); % energy of each sub image
en = en/sum(en(:)); % normalized by the total energy
q_ratio = 1./sqrt(en); % q ratios are inversely proportional to the square root of the energies
q_ratio = q_ratio./q_ratio(N,N) % normalising the q_ratios
q = quantratio1(Y,400, q_ratio, N, rise1);

z = quantratio2(q, q_ratio, N, 400, rise1);

dctbpp(regroup(z,N)/N, 16)

Z = lbt_dec(z, N, s);

draw(Z);


%% Investigating breaking into subimages
m = 2; %number of sub-images
w = size(X)/m;
X_sub = zeros([m, m, w]);
X_recon = zeros(256,256);
ens = zeros(m,m);
E = @(s) sum(s(:).^2);
log2 = @(s) log(s)/log(2);
N_level = log(w(1))/log(2) - 1;
for r = 1:m
    for c = 1:m
        X_sub(r,c,:,:) = X(1+(r-1)*w:r*w, 1+(c-1)*w:c*w);
        ens(r,c) = E(X_sub(r,c,:));
    end
end
ens = floor(40960*ens/sum(ens(:)));
for r = 1:m
    for c = 1:m
        X_s = squeeze(X_sub(r,c,:,:));
        [~,~,X_recon(1+(r-1)*w:r*w, 1+(c-1)*w:c*w),~] = lbt_opt_enc(X_s,4,16,1,sqrt(2),true,16,0.1,0,ens(r,c), m^2); %dwt_opt_enc(X_s,N_level,-1,0.5,-1,-1,1,16,ens(r,c));
    end
end
%%
figure(2)
draw(X_recon);
ssim(X-128, X_recon)
%[ssimval, rmsError, Z, q_opt] = lbt_opt_enc(X, N, M, rise1, s, opthuff, dcbits, ratio, ratio2, targ,m);

%% Just 3 level lbt
s = sqrt(2);
N = 4;
M = 16;
rise1=1;
[ssimval, rmsError, Z, q_opt] = lbt_opt_enc(X, N, M, rise1, s, true, 16, 0.1, 0, 40960,1);

