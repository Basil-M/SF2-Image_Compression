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
yticks(1:7);
yticklabels({'2','4','8','16','32','64','128'});
xlabel('Number of levels')
xticks(1:7);
zlabel('SSIM')

%best seems to always be M = 2^N
%in this case, 5|5?
%[s,n,xr,q,vlc] = dwt_opt_enc(X,5);
%draw(xr);


%% Different filters
dwtmode('per');
wname = 'bior3.3';
q0 = 1;
Y = dwt_f(X, 4, wname);

%C = Y{1};
%C = quantise(C,50);
%Y{1} = C;
%Yq2 = Y;

q_rf = dwt_q_ratios_f(Y, wname);
Yq1 = dwtquant1_f(Y, q0*q_rf, 0.5);
Yq2 = dwtquant2_f(Yq1, q0*q_rf, 0.5);

Xq = idwt_f(Yq2, wname);
ssim(Xq+128,X)

%% Compare high frequencies to the 
%generate lowpass filter
f_size = 10;
sig = 0.5;
f2 = ceil(f_size/2);
[X1,X2] = meshgrid(-f2:1:f2, -f2:1:f2);
F = mvnpdf([X1(:) X2(:)], [0 0], sig*[1 0; 0 1]);
F = reshape(F, 2*f2+1, 2*f2+1);

x_lp = conv2(X, F, 'same');
x_hp = X - x_lp;
[s_d,n,xr,q,vlc] = dwt_opt_enc(X,7);
evalc('[s_l,~,xrlbt,~] = jpegencdeclbtnlev(X, 4, 16, r, sqrt(2), true, 16, 0.1, 0);');
[~, ssimMap] = ssim(xr, X);
[s, ssimMapLBT] = ssim(xrlbt, X-128);

s_draw = ssimMap - min(ssimMap(:));
s_draw_l = ssimMapLBT - min(ssimMapLBT(:));

s_draw = s_draw/max(s_draw(:));
s_draw_l = s_draw_l/max(s_draw_l(:));

x_draw = x_hp - min(x_hp(:));

s_draw = s_draw*max(x_draw(:));
s_draw_l = s_draw_l*max(x_draw(:));

% draw(beside(x_draw, s_draw_l));
% figure(2);
% draw(beside(s_draw, x_draw));
draw(beside(beside(s_draw,x_draw), s_draw_l));
%% INVESTIGATING USING LBT

S_MAT = -Inf*ones([2,7]);
for n = 1:7
    [S_MAT(1,n),~,~,~,~] = dwt_opt_enc(X, n);
    [S_MAT(2,n),~,~,~,~] = dwt_opt_enc(X, n, -1, 0.5, 0, 1);
end

plot(S_MAT(1,:));
hold on;
plot(S_MAT(2,:));
ylabel('SSIM')
xlabel('Number of levels')
legend('No LBT','LBT on final DWT lowpass');

%% INVESTIGATING RISE
rises = 0.1:0.05:1.5;
S_MAT = -Inf*ones([2,length(rises)]);
%Q_MAT = -Inf*ones([2,length(rises)]);
i = 1;
for r = rises
   evalc('[S_MAT(1,i),~,~,~] = jpegencdeclbtnlev(X, 4, 16, r, sqrt(2), true, 16, 0.1, 0);');
   [S_MAT(2, i),~,~, ~, ~] = dwt_opt_enc(X,7, -1, r);
   i = i + 1;
end

plot(rises,S_MAT(1,:));
hold on;
plot(rises, S_MAT(2,:));
ylabel('SSIM value');
legend('LBT','DWT');
xlabel('Rise');
%% Experiment with resizing
scl = 1:0.2:2;
S_MAT = -Inf*ones([1, length(scl)]);
Q_MAT = S_MAT;
i = 1;
for sc = scl
    [~, ~,  X_rec, Q_MAT(i),~] = dwt_opt_enc(imresize(X, sc), 7);
    X_rec = imresize(X_rec, size(X,1)/size(X_rec,1));
    S_MAT(i) = ssim(X_rec, X);
    i = i + 1;
end
%failed

%% Experiment with SVD to remove unimportant detail
[U,S,V] = svd(X);
per_keep = SS(1,:);
S_MAT = -Inf*ones([length(per_keep) 1]);
i = 1; 
for p = per_keep
    S_p = S;
    ind = floor(p*length(S)/100);
    S_p(ind:length(S), ind:length(S)) = 0;
    X_p = U*S_p*V';
    [~,~,xr,~,~] = dwt_opt_enc(X_p,7);
    S_MAT(i) = ssim(xr, X);
    i = i + 1;
end


%% Investigating breaking into sub-images
m = 2; %number of sub-images
E = @(s) sum(s(:).^2);
log2 = @(s) log(s)/log(2);
S_MAT = [];
for m = [1 2]
    ens = zeros(m,m);
    X_recon = zeros(256,256);
    w = size(X,1)/m;
    X_sub = zeros([m, m, w, w]);
    N_level = log2(w) - 1;
    for r = 1:m
        for c = 1:m
            X_sub(r,c,:,:) = X(1+(r-1)*w:r*w, 1+(c-1)*w:c*w);
            ens(r,c) = E(X_sub(r,c,:));
        end
    end
    ens = floor(40960*ens/sum(ens(:)))
    for r = 1:m
        for c = 1:m
            X_s = squeeze(X_sub(r,c,:,:));
            [~,n,X_recon(1+(r-1)*w:r*w, 1+(c-1)*w:c*w),~,~] = dwt_opt_enc(X_s,N_level,-1,0.5,-1,-1,1,16,ens(r,c));
            %redistribute unused bit allowance to next one
            if ~(r == m && c == m)
                if c == m
                    ens(r + 1, c) = ens(r + 1, c) + (ens(r,c) - n);
                else
                    ens(r, c + 1) = ens(r, c + 1) + (ens(r,c) - n);
                end
            end    
        end
    end
    S_MAT = [S_MAT ssim(X_recon, X)];
end

%% Compress high and low freq separately
%generate lowpass filter
f_size = 10;
sig = 0.5;
f2 = ceil(f_size/2);
[X1,X2] = meshgrid(-f2:1:f2, -f2:1:f2);
F = mvnpdf([X1(:) X2(:)], [0 0], sig*[1 0; 0 1]);
F = reshape(F, 2*f2+1, 2*f2+1);

x_lp = conv2(X, F, 'same');
x_hp = X - x_lp;

BW1 = edge(X, 'sobel');
dil = imdilate(BW1, strel('disk',4,4));

X_edges = X;
X_edges(~dil) = 0;
X_noedges = X;
X_noedges(dil) = 0;

bits_edges = 10000;
[~,bits_edges,xr_edge,~,~] = dwt_opt_enc(X_edges,7,-1,0.5,-1,-1,1,16,bits_edges);
[~,~,xr_noedge,~,~] = dwt_opt_enc(X_noedges,7,-1,0.5,-1,-1,1,16,40960 - bits_edges);

xr = xr_noedge + xr_edge;
