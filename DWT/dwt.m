function Y = dwt(X,h1,h2)

%DWT Discrete Wavelet Transform
%  Y = DWT(X, H1, H2) returns a 1-level 2-D discrete wavelet
%  transform of X.
%
%  If filters H1 and H2 are given, then they are used,
%  otherwise the LeGall filter pair are used.

if nargin < 3,
  h1=[-1 2 6 2 -1]/8;
  h2=[-1 2 -1]/4;
end

[m,n] = size(X);
Y = zeros(m,n);

n2 = n/2;
t = 1:n2;
Y(:,t) = rowdec(X,h1);
Y(:,t+n2) = rowdec2(X,h2);

X = Y';
m2 = m/2;
t = 1:m2;
Y(t,:) = rowdec(X,h1)';
Y(t+m2,:) = rowdec2(X,h2)';



%{
per_keep = 10;
[U,S,V] = svd(Y);
ind = ceil(per_keep*m/100);
S(ind:m, ind:m) = 0;
Y = U*S*V';

for r = 1:2
    for c = 1:2
        %don't remove detail from lowpass as this will propagate through
        %subsequent layers
        if ~(r==1 && c == 1)
            Y_s = Y(1+(r-1)*m2:r*m2, 1+(c-1)*m2:c*m2);   
            [U,S,V] = svd(Y_s);
            ind = ceil(per_keep*m2/100);
            S(ind:length(S), ind:length(S)) = 0;
            Y(1+(r-1)*m2:r*m2, 1+(c-1)*m2:c*m2) = U*S*V';
        end
    end
end
%}
end

