function [Y] = dwt_f(X,N,wname)
%DWT_F Summary of this function goes here
%   Detailed explanation goes here

if exist('wname','var')
    if wname~=-1
        Y = cell(1,2);
        [Y{1},Y{2}] = wavedec2(X, N, wname);
    else
        Y = idwt(X);
    end
else
    Y = idwt(X);
end

end

