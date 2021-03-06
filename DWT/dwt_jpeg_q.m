function [q_opt] = dwt_jpeg_q(X,N_LEVELS, M,rise,N_sup,N_LBT, opthuff, targ)
%DWT_JPEG_Q Summary of this function goes here
%   Detailed explanation goes here

%set default parameter values
if~exist('rise','var'); rise = 0.5; end
if ~exist('opthuff','var'); opthuff = true; end
if ~exist('M','var'); M = -1; end
if ~exist('N_sup','var'); N_sup = -1; end
%if ~exist('wname','var'); wname = -1; end
if ~exist('N_LBT','var'); N_LBT = -1; end
vlc_q = @(q) dwt_enc(X, N_LEVELS, M, q, rise, N_sup, N_LBT, opthuff);
nbits = @(q) jpegbits(vlc_q(q), opthuff, true);

%fn = @(q) (nbits(q) - 40960)^2; %function to minimise
q_cur = 20;
precision = 4;
if ~exist('targ','var'); targ = 40960; end
step = 1;
%initialisation
while(nbits(q_cur) > targ) && (q_cur < 256)
    q_cur = q_cur + 10;
end

if(q_cur == 260)
    disp('Could not optimise q for these parameters - optimal q-step too large.\n');
    q_opt = -1;
else
    for p = 1:precision
       while nbits(q_cur) < targ
        q_cur = q_cur - step;
        if(q_cur <= 0)
            q_opt = -1;
            return 
        end
       end
       q_cur = q_cur + step;
       step = step/10; 
    end

    q_opt = q_cur + step;
end
end


