% doing the same as jpegencdeclbt but using the energies of the sub images
% to work out how much to encode each image
function [ssimval, rmsError, Z] = jpegencdeclbtenergy(X, N, M, rise1, s, opthuff, dcbits);

%find the step size which gives 5kB
q_opt = lbtjpegqenergy(X, N, M, rise1, s);

% encode the image
%[vlc bits huffval] = jpegenclbt(X, q_opt, N, M, rise1, s, opthuff, dcbits);
[vlc, bits, huffval, q_ratio] = jpegenclbtenergy(X, q_opt, N, M, rise1, s, opthuff, dcbits);


%decode the image
%Z = jpegdeclbt(vlc, q_opt, N, M, rise1, s, bits, huffval, dcbits);
Z = jpegdeclbtenergy(vlc, q_opt, q_ratio, N, M, rise1, s, bits, huffval, dcbits);

% compare the decoded image to the original
ssimval = ssim(Z,X-128); 

% rms error
rmsError = std(X(:)-Z(:));
end 

% finds the value of q_opt which makes the jpeg encoded image 5 Kbs
function q_opt = lbtjpegqenergy(X, N, M, rise1, s);

target = 40960; % Max number of bits
f = @(q) jpegbitsenergy(jpegenclbtenergy(X, q, N, M, rise1, s, true, 64),true, false, N);% - target).^2;
q_opt = 500;
while f(q_opt) > target
    q_opt = q_opt + 5;
end



precision = 4;
i = 0;
step = 1;
while i<precision
    q_opt = q_opt - step;
    if f(q_opt) > target %then want to reduce f by increasing q
        q_opt = q_opt + step;
        step = step/10;
        i = i+1;
    end
end


end

function nbits = jpegbitsenergy(vlc, opthuff, dwt, N)
%JPEGBITS Calculates number of bits of a huffman code
%   Detailed explanation goes here
if opthuff
    nbits = vlctest(vlc) + 1424;
else
    nbits = vlctest(vlc);
end

if dwt
    nbits = nbits + 88;
else
    nbits = nbits + 120 + N^2*32;
end

end



