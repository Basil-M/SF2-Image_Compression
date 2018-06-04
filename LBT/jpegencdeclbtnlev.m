function [ssimval, rmsError, Z, q_opt] = jpegencdeclbtnlev(X, N, M, rise1, s, opthuff, dcbits, ratio, ratio2);

%find the step size which gives 5kB
q_opt = lbtjpegqnlev(X, N, M, rise1, s, opthuff, dcbits, ratio, ratio2);

% encode the image
[vlc, bits, huffval,~,~] = jpegenclbtnlev(X, q_opt, N, M, rise1, s, opthuff, dcbits, ratio,ratio2);

%decode the image
Z = jpegdeclbtnlev(vlc, q_opt, N, M, rise1, s, bits, huffval, dcbits, ratio,ratio2);

% compare the decoded image to the original
ssimval = ssim(Z,X); 

% rms error
rmsError = std(X(:)-Z(:));
end

% finds the value of q_opt which makes the jpeg encoded image 5 Kbs
function q_opt = lbtjpegqnlev(X, N, M, rise1, s, opthuff, dcbits, ratio, ratio2);

target = 40960; % Max number of bits
f = @(q) jpegbits(jpegenclbtnlev(X, q, N, M, rise1, s, opthuff, dcbits, ratio, ratio2),true, false);% - target).^2;
q_opt = 40;
while f(q_opt) > target
    q_opt = q_opt + 5;
end



precision = 5;
i = 0;
step = 1;
while i<precision
    q_opt = q_opt - step;
    fopt = f(q_opt);
    if f(q_opt) > target %then want to reduce f by increasing q
        q_opt = q_opt + step;
        step = step/10;
        i = i+1;
    end
end


end
