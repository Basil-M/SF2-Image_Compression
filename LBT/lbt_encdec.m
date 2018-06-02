function [rmsError, totalBits, compRatio, Zp] = lbt_encdec(X, N, q, s);
I = size(X, 1); % I am assuming size is just the size of  one edge. Check this
[Pf, Pr] = pot_ii(N,s); 

t = (1+N/2):(I-N/2);% N is the DCT size, I is the image size
Xp = X; % Do I need to give X zero mean?
Xp(t,:) = colxfm(Xp(t,:),Pf);
Xp(:,t) = colxfm(Xp(:,t)',Pf)';

%[~,totalBits,compRatio,Z] = dct_encdec(Xp, N, q);

CN = dct_ii(N);
Xpzero = Xp-128;
Y=colxfm(colxfm(Xpzero,CN)',CN)';
% quantise Y
Yq = quantise(Y, q);
% reconstruct the image from Yq
Zp = colxfm(colxfm(Yq',CN')',CN');

Yr = regroup(Yq, N);

%Use below for LBT
totalBits = dctbpp(Yr, 16);
Xq = quantise(X, 17);
refBits = bpp(Xq)*numel(Xq);
compRatio = refBits/totalBits;


Zp(:,t) = colxfm(Zp(:,t)',Pr')';
Zp(t,:) = colxfm(Zp(t,:), Pr');
rmsError = std(X(:)-Zp(:));


