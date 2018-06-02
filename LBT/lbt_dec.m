%decode the lbt.
function Zp = lbt_dec(Y, N, s);
%size of the image
I = size(Y, 1);

t = (1+N/2):(I-N/2);% N is the DCT size, I is the image size

CN = dct_ii(N);
% get the post filter
[~, Pr] = pot_ii(N,s); 
% reconstruct the image from Yq
Zp = colxfm(colxfm(Y',CN')',CN');
% post filter the image
Zp(:,t) = colxfm(Zp(:,t)',Pr')';
Zp(t,:) = colxfm(Zp(t,:), Pr');