function y = quantratio2(q, q_ratio, N, step, rise1)

% QUANT2 Reconstruct matrix from quantised values
%  Y = QUANT2(Q, step, rise1) Reconstructs the matrix Y from integers
%  Q using steps of width step.
%  q_ratio is the ratio in which each sub image should be quantised
%  N is the size of the dct transform
%
%  The result is the reconstructed values. If rise1 is defined, the first
%  step rises at rise1, otherwise it rises at step/2 to give a uniform
%  quantiser with a step centred on zero.
%  In any case the quantiser is symmetrical about zero.

if step <= 0, y = q; return, end

if nargin <= 4, rise1 = step/2; end

% Reconstruct quantised values and incorporate sign(q).
%y = q * step + sign(q) * (rise1 - step/2);
for i=1:4
    for j=1:4
        qsub = q(i:N:256,j:N:256);
        y(i:N:256,j:N:256) = quant2(qsub, step*q_ratio(i,j), rise1);%qsub * step * q_ratio(i,j) + sign(qsub) * (rise1 - step*q_ratio(i,j)/2);
    end
end
