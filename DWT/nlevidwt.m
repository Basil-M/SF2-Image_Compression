%multi level inverse DWT
function Xr = nlevidwt(Y, N);
max_size = size(Y,1);
m=max_size/(2.^N);
Xr = Y;
while m<max_size
    m = m*2;
    t = 1:m;
    % reconstruct each layer of sub-images. Ending by reconstructing
    % 256x256 image. 
    Xr(t,t) = idwt(Xr(t,t));
end