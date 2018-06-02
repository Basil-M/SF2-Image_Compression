% finds the value of q_opt which makes the jpeg encoded image 5 Kbs
function q_opt = lbtjpegq(X, N, M, rise1, s);

target = 40960; % Max number of bits
f = @(q) jpegbits(jpegenclbt(X, q, N, M, rise1, s, true),true, false);% - target).^2;
q_opt = 30;
while f(q_opt) > target
    q_opt = q_opt + 5;
end



precision = 4;
i = 0;
step = 1;
while i<precision
    q_opt = q_opt - step;
    if f(q_opt) > target %then want to reduce f by increasing q
        q_opt = q_opt + step;
        step = step/10;
        i = i+1;
    end
end


end


