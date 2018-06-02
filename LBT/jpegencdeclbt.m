function [ssimval, rmsError, Z] = jpegencdeclbt(X, N, M, rise1, s, opthuff, dcbits);

%find the step size which gives 5kB
q_opt = lbtjpegq(X, N, M, rise1, s);

% encode the image
[vlc bits huffval] = jpegenclbt(X, q_opt, N, M, rise1, s, opthuff, dcbits);

%decode the image
Z = jpegdeclbt(vlc, q_opt, N, M, rise1, s, bits, huffval, dcbits);

% compare the decoded image to the original
ssimval = ssim(Z,X-128); 

% rms error
rmsError = std(X(:)-Z(:));