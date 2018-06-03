%encodes with lbt and supressed the highest frequency components
function Y = lbt_encsup(X, N, s);
%size of the image
I = size(X, 1);
% pre filter
[Pf, ~] = pot_ii(N,s); 

t = (1+N/2):(I-N/2);% N is the DCT size, I is the image size
Xp = X; 
Xp(t,:) = colxfm(Xp(t,:),Pf);
Xp(:,t) = colxfm(Xp(:,t)',Pf)';

CN = dct_ii(N);
%supress the final row and column of sub images
CN(N,:) = 0;
%Giving Xp zero mean
Xpzero = Xp-128;
Y=colxfm(colxfm(Xpzero,CN)',CN)';

