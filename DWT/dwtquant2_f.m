function Yq = dwtquant2_f(Yq1,dwtstep, rise)

%DWTQUANT1 Quantise DWT and calculate entropy/number of bits
%   Y = encoded DWT
%   N = Number of levels
%   dwtstep = Array of quantisation steps

%  The result is the quantised integers Q. If rise1 is defined,
%  the first step rises at rise1, otherwise it rises at step/2 to
%  give a uniform quantiser with a step centred on zero.
%  In any case the quantiser is symmetrical about zero.

%Make fake Y
S = Yq1{2};
C = Yq1{1};
N = length(S) - 2;
Yq = Yq1;
%start with approximation matrix
q_ind = 1;          %current quantisation step
cur_ind = 1;        %stores current starting index
len_n = prod(S(1,:)) - 1;
C(cur_ind:cur_ind+len_n) = quant1(C(cur_ind:cur_ind+len_n),dwtstep(1),rise*dwtstep(1));
q_ind = q_ind + 1;
cur_ind = cur_ind + len_n + 1;
%inverse order to original dwt_q_ratios function
for i = 1:N   
    len_n = prod(S(i+1,:)) - 1;     %number of elements in C column vector this corresponds to
    for k = 1:3
        fprintf("Unquantising range %i to %i with step %0.2f\n", cur_ind, cur_ind+len_n, dwtstep(q_ind));
        C(cur_ind:cur_ind+len_n) = quant2(C(cur_ind:cur_ind+len_n),dwtstep(q_ind),rise*dwtstep(q_ind));
        q_ind = q_ind + 1;
        cur_ind = cur_ind + len_n + 1;
    end
end
Yq{1} = C;
end
