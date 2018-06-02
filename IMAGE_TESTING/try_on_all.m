function [output] = try_on_all(f, N)
%TRY_ON_ALL Try function f on all images
%   Optional: N, number of images to try it on

fol_inf = dir('Images');
fol_inf(~[fol_inf.isdir])=[];
output = [];
m = 1;

for i = 3:length(fol_inf)
    p = dir(['Images/' fol_inf(i).name]);
    if nargin == 2
        k_max = min(3+ ceil(N/(length(fol_inf)-2)), length(p));
    else
        k_max = length(p);
    end    
    for k = 3:k_max
        f_name = ['Images/' fol_inf(i).name '/' p(k).name];
        output(m,:) = f(image_reader(f_name));
        m = m + 1;
    end
end

end

