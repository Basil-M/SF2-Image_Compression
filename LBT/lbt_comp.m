%Identical to lbt_encdec.m but just returns the compression ratio so this
%can be optimized wrt s.
function [compRatio, q] = lbt_comp(X, N, s);
I = size(X, 1); % I am assuming size is just the size of  one edge. Check this
[Pf, Pr] = pot_ii(N,s); 

t = (1+N/2):(I-N/2);% N is the DCT size, I is the image size
Xp = X; % Do I need to give X zero mean?
Xp(t,:) = colxfm(Xp(t,:),Pf);
Xp(:,t) = colxfm(Xp(:,t)',Pf)';
q = opStep(X, N, 3, s);

[~, ~, compRatio,~] = dct_encdec(Xp, N, q);

%{
Zp(:,t) = colxfm(Zp(:,t)',Pr')';
Zp(t,:) = colxfm(Zp(t,:), Pr');
%}

