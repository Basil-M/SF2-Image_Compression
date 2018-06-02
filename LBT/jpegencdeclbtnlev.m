function [ssimval, rmsError, Z, q_opt] = jpegencdeclbtnlev(X, N, M, rise1, s, opthuff, dcbits);

%find the step size which gives 5kB
q_opt = lbtjpegqnlev(X, N, M, rise1, s, opthuff, dcbits);

% encode the image
[vlc, bits, huffval] = jpegenclbtnlev(X, q_opt, N, M, rise1, s, opthuff, dcbits);

%decode the image
Z = jpegdeclbtnlev(vlc, q_opt, N, M, rise1, s, bits, huffval, dcbits);

% compare the decoded image to the original
ssimval = ssim(Z,X-128); 

% rms error
rmsError = std(X(:)-Z(:));
end

% finds the value of q_opt which makes the jpeg encoded image 5 Kbs
function q_opt = lbtjpegqnlev(X, N, M, rise1, s, opthuff, dcbits);

target = 40960; % Max number of bits
f = @(q) jpegbits(jpegenclbtnlev(X, q, N, M, rise1, s, opthuff, dcbits),true, false);% - target).^2;
q_opt = 30;
while f(q_opt) > target
    q_opt = q_opt + 5;
end



precision = 5;
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
