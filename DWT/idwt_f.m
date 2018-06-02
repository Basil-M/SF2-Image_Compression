function X = idwt_f(Y, wname)
%IDWT_F Summary of this function goes here
%   Detailed explanation goes here

if exist('wname','var')
    if wname~=-1
        X = waverec2(Y{1}, Y{2}, wname);
    else
        Y = idwt(X);
    end
else
    Y = idwt(X);
end

end

