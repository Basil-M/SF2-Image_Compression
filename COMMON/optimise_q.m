function [q_opt] = optimise_q(fn, targ)
%OPTIMISE_Q Summary of this function goes here
%   Detailed explanation goes here

f = @(q) (fn(q) - targ).^2;

% Find best staring point via grid search
c_min = Inf;
for q = 5:1:50
    f_val = f(q);
    if f_val < 0.01^2
        q_opt = q;
        break;
    elseif f_val < c_min
        c_min = f_val;
        q_opt = q;
    end
end

[q_opt, ~] = fminsearch(f,q_opt);% , options);

end

