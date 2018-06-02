%encodes with lbt and NO quantisation
function Y = lbt_justenc(X, N, s);
%size of the image
I = size(X, 1);
% pre filter
[Pf, ~] = pot_ii(N,s); 

t = (1+N/2):(I-N/2);% N is the DCT size, I is the image size
Xp = X; 
Xp(t,:) = colxfm(Xp(t,:),Pf);
Xp(:,t) = colxfm(Xp(:,t)',Pf)';

CN = dct_ii(N);
%Giving Xp zero mean % not doing this anymore. Just zero mean original
%image
%Xpzero = Xp;
Y=colxfm(colxfm(Xp,CN)',CN)';