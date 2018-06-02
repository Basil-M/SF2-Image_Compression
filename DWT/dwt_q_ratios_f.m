function [dwtstep] = dwt_q_ratios_f(Y, wname)
%DWT_Q_RATIOS Finds equal MSE Q ratios for n level dwt
%   sz = Size of image
%   N = number of levels of dwt ratio

%Make fake Y
S = Y{2};
C = Y{1};
N = length(S) - 2;
C_zeros = zeros(size(C));
%Set impulse value
E_imp = 100;

E = @(Z) sum(Z(:).^2);  %Calculate energy
E_rats = zeros(3*N + 1,1);

%start with approximation matrix
e_ind = 1;          %stores index of energy matrix
C = C_zeros;
cur_ind = 1;        %stores current starting index
len_n = prod(S(1,:));
C(cur_ind + ceil(len_n/2)) = E_imp;

%reconstruct
Y{1} = C;
X_t = idwt_f(Y, wname);
E_rats(1) = E(X_t)/E_imp^2;
e_ind = e_ind + 1;
cur_ind = cur_ind + len_n;
%inverse order to original dwt_q_ratios function
for i = 1:N   
    len_n = prod(S(i+1,:));     %number of elements in C column vector this corresponds to
    l2 = ceil(len_n/2);
    for k = 1:3
        C = C_zeros;
        C(cur_ind+l2) = E_imp;
        Y{1} = C;
        X_t = idwt_f(Y, wname);
        E_rats(e_ind) = E(X_t)/E_imp^2;
        e_ind = e_ind + 1;
        cur_ind = cur_ind + len_n;
    end
end

dwtstep = 1./sqrt(E_rats);
dwtstep = dwtstep/max(dwtstep); %normalise
end

