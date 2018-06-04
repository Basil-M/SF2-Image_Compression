% DWT with N layers 
function Y =nlevdwt(X, N)
m = size(X, 1); 
% not sure this is correct place to subtract 128 but gives good image
Y = dwt(X);%dwt(X-128);
i = 1;
while i<N
    m=m/2;
    t=1:m;
    Y(t,t) = dwt(Y(t,t));
    i = i+1;
end

end
