function [ssim_r, nbits, Xr, Q, vlc] = dwt_opt_enc(X, N_LEVELS, M, rise, N_sup, N_LBT, opthuff, dcbits, targ)
%DWT_OPT_ENC Optimises Q and returns encoded dwt in one (for convenience!)
%dwt_opt_enc(X, 7,-1,-1,-1,'haar')
if ~exist('rise','var')
    rise = 0.5;
elseif rise < 0
    rise = 0.5;
end

if ~exist('N_sup', 'var')    N_sup = 0;

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
    M = -1;
end

if ~exist('dcbits','var')
    dcbits =16;
end

if ~exist('N_LBT','var')
    N_LBT = 0;
end

fprintf('Optimising q for N = %i, M = %i, rise = %0.2f, N_sup=%i\n',N_LEVELS, M, rise, N_sup);
if exist('targ','var')
    Q = dwt_jpeg_q(X, N_LEVELS, M, rise, N_sup,N_LBT, opthuff, targ);
else
    Q = dwt_jpeg_q(X, N_LEVELS, M, rise, N_sup,N_LBT, opthuff);
end
if(Q <= 0)
    Xr = [];
    nbits = [];
    ssim_r = -Inf;
    vlc = [];
    bits = [];
    huffval = [];
else
    fprintf('Final encoding with optimal q0 = %0.2f\n',Q);
    [vlc, bits, huffval] = dwt_enc(X, N_LEVELS, M, Q, rise, N_sup,N_LBT, opthuff, dcbits);

    Xr = dwt_dec(vlc, N_LEVELS, M, Q, rise,N_LBT, bits, huffval, dcbits, size(X,1), size(X,2));
    %Xr = dwt_dec(vlc, N_LEVELS, M, Q*dwt_q_ratios_X(X, N_LEVELS), rise,N_LBT, bits, huffval, dcbits, size(X,1), size(X,2));

    nbits = jpegbits(vlc, opthuff, true);
    
    ssim_r = max([ssim(Xr, X), ssim(Xr, X - 128), ssim(Xr, X + 128)]);
    
end
end

