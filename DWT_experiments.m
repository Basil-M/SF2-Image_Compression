N_DWT = 7;

inBits = @(kB) kB*8*1024;
inkB = @(Bits) Bits/(8*1024);

%% BASIC DWT
% This is the basic DWT
% By default it uses a rise of 1, no suppression, subimage width of 2^N_DWT
[s,n,xr,q,vlc] = dwt_opt_enc(X,N_DWT);

fprintf("Image takes up %0.2fkB and has an SSIM of %0.4f", inkB(n), s);
figure('Name','Reconstructed image'); draw(xr);

%% INVESTIGATING SUPPRESSION
recon_X = zeros([3*N_DWT + 1, size(X)]);
S = []; nbits = []; Q = [];
for n = 0:(3*N_DWT - 1)
    [S(n + 1),nbits(n+1),recon_X(n+1, :, :),Q(n+1),~] = dwt_opt_enc(X, N_DWT, 2^N_DWT, 1, n);
end

plot(0:(3*N_DWT - 1),S,'o');
xlabel('Number of suppressed sub-images');
ylabel('SSIM');
yyaxis right;
plot(0:(3*N_DWT - 1),Q,'x');
ylabel('Quantisation q0');

% every third image
Xdraw = beside(X,squeeze(recon_X(1,:,:)));
for n = 4:3:(3*N_DWT -  1)
    Xdraw = beside(Xdraw, squeeze(recon_X(n,:,:)));
end
figure('Name','Suppressed reconstructions');
draw(Xdraw);

%% INVESTIGATING SUB-IMAGE WIDTH

S_MAT = -Inf*ones([7,7]);
for n = 1:7
    for m = n:7
        [S_MAT(n, m),~,~,~,~] = dwt_opt_enc(X, n, 2^m);
    end
end

surf(S_MAT);
ylabel('Encoding width')
yticklabels({'1','2','4','8','16','32','64','128','256'});
xlabel('Number of levels')
xticks(1:7);
zlabel('SSIM')