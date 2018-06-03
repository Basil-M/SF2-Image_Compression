function q = quantratio1(x, step, q_ratio, N, rise1)

% QUANT1 Quantise a matrix
%  Q = QUANT1(X, step, q_ratio, rise1) quantises the matrix X using steps
%  of width step*q_ratio.
%  N is the size of the dct transform
%  q_ratio is the proportion of the step with which to quantise each of the subimages with 
%  The result is the quantised integers Q. If rise1 is defined,
%  the first step rises at rise1, otherwise it rises at step/2 to
%  give a uniform quantiser with a step centred on zero.
%  In any case the quantiser is symmetrical about zero.
% quantises each sub image of the lbt with a different quantisation step
% according to q_ratio.

if step <= 0, q = x; return, end

if nargin <= 4, rise1 = step/2; end

% Quantise abs(x) to integer values, and incorporate sign(x)..
%q = max(0,ceil((abs(x) - rise1)/step)) .* sign(x);
for i=1:4
    for j=1:4
        xsub = x(i:N:256,j:N:256);
        q(i:N:256,j:N:256) =  quant1(xsub, step*q_ratio(i,j), rise1);% max(0, ceil(abs(xsub) - rise1)/(step*q_ratio(i,j))) .* sign(xsub);
    end
end
