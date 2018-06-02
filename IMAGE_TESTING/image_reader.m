function [X] = image_reader(f_path)
%IMAGE_READER Reads an image X, converts to 256 x 256 grayscale
%   Detailed explanation goes here
X = imread(f_path);
if ndims(X) ~= 2
    X = rgb2gray(X);
end

%% PAD WITH ZEROS
% X = imresize(X, 256/max(size(X)));
% 
% 
% if mod(min(size(X)),2) ~= 0
%     pad_amount = (256 - (min(size(X)) - 1))/2; %pad extra row/col
%     odd_pad = true;                          %note to delete one row/col later
% else
%     pad_amount = (256 - min(size(X)))/2;
%     odd_pad = false;
% end
% 
% if size(X,1) > size(X,2)
%     %more rows than columns
%     X = padarray(X, [0 pad_amount]);
% else
%     X = padarray(X, [pad_amount 0]);
% end
% if odd_pad
%     X = X(1:256, 1:256);
% end

%% RESIZE AND CROP
X = imresize(X, 256/min(size(X)));
if max(size(X)) ~= 256
    s_max = ceil((max(size(X)) - 256)/2);
    if size(X, 1) > size(X, 2)
        %more rows than columns
        X = X(s_max:s_max + 255, :);
    else
        X = X(:, s_max:s_max + 255);
    end
end
size(X);
X = double(X);


end

