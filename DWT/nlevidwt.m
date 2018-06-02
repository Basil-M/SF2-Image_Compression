%multi level inverse DWT
function Xr = nlevidwt(Y, N);

m=256/(2.^N);
Xr = Y;
while m<256
    m = m*2;
    t = 1:m;
    % reconstruct each layer of sub-images. Ending by reconstructing
    % 256x256 image. 
    Xr(t,t) = idwt(Xr(t,t));
end