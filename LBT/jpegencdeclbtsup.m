%  supressing the high frequency components in the encoding
function [ssimval, rmsError, Z] = jpegencdeclbtsup(X, N, M, rise1, s, opthuff, dcbits);

%find the step size which gives 5kB
q_opt = lbtjpegq(X, N, M, rise1, s);

% encode the image
[vlc bits huffval] = jpegencsup(X, q_opt, N, M, rise1, s, opthuff, dcbits);

%decode the image
Z = jpegdeclbt(vlc, q_opt, N, M, rise1, s, bits, huffval, dcbits);

% compare the decoded image to the original
ssimval = ssim(Z,X-128); 

% rms error
rmsError = std(X(:)-Z(:));
end

% finds the value of q_opt which makes the jpeg encoded image 5 Kbs
function q_opt = lbtjpegq(X, N, M, rise1, s);

target = 40960; % Max number of bits
f = @(q) jpegbits(jpegencsup(X, q, N, M, rise1, s, true),true, false);% - target).^2;
q_opt = 30;
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


