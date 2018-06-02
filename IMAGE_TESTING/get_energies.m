%[ssimval, rmsError, Z, q_opt] = jpegencdeclbtnlev(X, 4, 15, 1, sqrt(2), true, 16, 0.2, 0);

N_DCT = 16;
C = dct_ii(8);

Y = @(X) lbt_enc(X, N_DCT, 1, sqrt(2), 1);
en = @(X) dct_energies(Y(X), N_DCT);

N_try = 1;

fol_inf = dir('IMAGE_TESTING/IMAGES');
fol_inf(~[fol_inf.isdir])=[];
output = [];
m = 1;
f_name = cell(0);
for i = 3:length(fol_inf)
    p = dir(['IMAGE_TESTING/IMAGES/' fol_inf(i).name]);
    if exist('N_try','var')
        k_max = min(3+ ceil(N_try/(length(fol_inf)-2)), length(p));
    else
        k_max = length(p);
    end    
    for k = 3:k_max
        f_name{m} = ['IMAGE_TESTING/IMAGES/' fol_inf(i).name '/' p(k).name];
        X_l = image_reader(f_name{m});
        [ssim(m, 2),~,~, q_opt(m,2), ~] = dwt_opt_enc(X, 7, -1, 0.7);
        [ssim(m, 1),~,~,q_opt(m,1)] = jpegencdeclbtnlev(X_l, 4, 15, 1, sqrt(2), true, 16, 0.2, 0);
        energies(m,:, :) =  en(X_l);
        m = m + 1;
    end
end

%% Get files from output
N_STORE = 10;
files = cell(N_DCT, N_DCT, 3, N_STORE);
for r = 1:N_DCT
    for c = 1:N_DCT
        [A, I] = sort(output(:, r, c));
        m_ind = floor(0.5*(length(I) - N_STORE));
        
        for n = 1:N_STORE
            %highest energy
            files{r, c, 1, n} = f_name{I(1 + n)};
            
            %medium energy
            files{r, c, 2, n} = f_name{I(m_ind + n)};
            
            %high energy
            files{r,c,3,n} = f_name{I(length(I) - N_STORE + n)};
        end
    end
end
 
        