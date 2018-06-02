%encodes with lbt. rise1 is the only optional parameter. The default value
%if 0.5. rise1 is the proportion of q that you want to use for step size
%around zero. It is multiplied by q before being passed into quantise.
function Yq = lbt_enc(X, N, q, s, rise1);
%size of the image
I = size(X, 1);
% pre filter
[Pf, ~] = pot_ii(N,s); 

t = (1+N/2):(I-N/2);% N is the DCT size, I is the image size
Xp = X; 
Xp(t,:) = colxfm(Xp(t,:),Pf);
Xp(:,t) = colxfm(Xp(:,t)',Pf)';

CN = dct_ii(N);
%Giving Xp zero mean
Xpzero = Xp-128;
Y=colxfm(colxfm(Xpzero,CN)',CN)';

%default value of rise1
if ~exist('rise1', 'var')
    rise1=0.5;
end
%convert proportion of step size to actual size of rise
rise1 = rise1*q;
% quantise Y
Yq = quantise(Y, q, rise1);
