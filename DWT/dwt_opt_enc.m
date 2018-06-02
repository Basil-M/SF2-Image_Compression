function [ssim_r, nbits, Xr, Q, vlc] = dwt_opt_enc(X, N_LEVELS, M, rise, N_sup, wname, opthuff, dcbits)
%DWT_OPT_ENC Optimises Q and returns encoded dwt in one (for convenience!)
%dwt_opt_enc(X, 7,-1,-1,-1,'haar')
if ~exist('rise','var')
    rise = 1;
elseif rise < 0
    rise = 1;
end

if ~exist('N_sup', 'var')
    N_sup = 0;
elseif N_sup < 0
    N_sup = 0;
end
if ~exist('wname','var')
    wname = -1;
end

if ~exist('opthuff','var')
    opthuff = true;
elseif opthuff == -1
    opthuff = true;
end

if ~exist('M','var')
    M = 2^N_LEVELS;
elseif M < 0
    M = 2^N_LEVELS;
end

if ~exist('dcbits','var')
    dcbits =16;
end

Q = dwt_jpeg_q(X, N_LEVELS, M, rise, N_sup, opthuff);

[vlc, bits, huffval] = dwt_enc(X, N_LEVELS, M, Q, rise, N_sup,opthuff, dcbits);

Xr = dwt_dec(vlc, N_LEVELS, M, Q, rise, bits, huffval, dcbits, size(X,1), size(X,2));

nbits = jpegbits(vlc, opthuff, true);

ssim_r = ssim(Xr, X);
end

