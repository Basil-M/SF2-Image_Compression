%% DWT encoding
% only optimisation is optimisation of the rise parameter

rises = [0.45 0.5 1];
rise = 0.5;
s_dwt = -1;
%{
q_opt = 0;
for r = rises
    evalc('[s,n,~,q,~] = dwt_opt_enc(X,7,-1,r);');
    if s > s_dwt
        q_opt = q;
        s_dwt = s;
        rise = r;
    end
end
%}

%% LBT encoding
% only optimisation is optimisation of the ratio parameter and the qstep

s_lbt = -1;
i=1;
ssimval =0;
q_opt_lbt_vec = 0;
for ratio1=0.05:0.05:0.35
    evalc('[ssimval(i), ~, ~, q_opt_lbt_vec(i)] = jpegencdeclbtnlev(X, 4, 16, 0.5, sqrt(2), true, 16, ratio1, 0);');
    i = i+1;
end

ratio1 = 0.05:0.05:0.35;

[s_lbt,I] = max(ssimval);

q_opt_lbt = q_opt_lbt_vec(I);
ratio = ratio1(I);


    
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
        
        [vlc, bits, huffval, q_opt_lbt, ratio] = jpegenclbtnlev(X, q_opt_lbt, 4, 16, 0.5, sqrt(2), true, 16, ratio, 0);

        
        
else
    %% Decoder needs variables
        % q_opt
        % rise
        % vlc
        % bits
        % huffval
    [vlc, bits, huffval] = dwt_enc(X, 7, -1, q_opt, rise);
end

%save(fname+'cmp');
save(strcat(fname, 'cmp'));
