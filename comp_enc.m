%% ADD PATHS
addpath('DWT');
addpath('COMMON');
addpath('LBT');

%% DWT encoding
% only optimisation is optimisation of the rise parameter
rises = 0.4:0.5:1.3;
rise = 0.5;
s_dwt = -1;
q_opt = 0;
for r = rises
    evalc('[s,n,~,q,~] = dwt_opt_enc(X,7,-1,r);');
    if s > s_dwt
        q_opt = q;
        s_dwt = s;
        rise = r;
    end
end

%% LBT encoding
% only optimisation is optimisation of the ratio parameter



s_lbt = -1;
    
%% CHOOSING BEST
if s_lbt > s_dwt
    % ENCODE USING LBT
    clear rise; % this deletes variable rise
                % decoder can tell using exist('rise','var') whether we
                % used DWT or LBT
                
    %% Encoder needs variables
        % q_opt
        % vlc
        % bits
        % huffval
        % and...?
else
    %% Decoder needs variables
        % q_opt
        % rise
        % vlc
        % bits
        % huffval
    [vlc, bits, huffval] = dwt_enc(X, 7, -1, q_opt, rise);
end

save(fname+"cmp");

