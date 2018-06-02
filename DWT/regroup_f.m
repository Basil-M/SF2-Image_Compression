function [Ydraw] = regroup_f(Y, wname)
%REGROUP_F Given bookkeeping matrix S and 
%   Detailed explanation goes here
C = Y{1}; S = Y{2};

below = @(x1, x2) rot90(beside(rot90(x1), rot90(x2)), 3);

%structure: A(N), H(N), V(N), D(N), H(N-1), V(N-1), D(N-1)
%S(1,:) = size of approximation coefficients
Ydraw = imresize(appcoef2(C,S, wname, 1),S(:,1));
N = length(S) - 2;
for i = 1:N
    [Hi, Vi, Di] = detcoef2('all',C,S,i);
    Ybot = beside(Hi, Di);
    Ydraw = below(beside(Ydraw, Vi),Ybot);
end
end

