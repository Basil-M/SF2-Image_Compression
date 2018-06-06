
%[ssimval, rmsError, Z, q_opt] = jpegencdeclbtnlev(X, 4, 15, 1, sqrt(2), true, 16, 0.2, 0);

N_DCT = 16;;

Y = @(X) lbt_enc(X, N_DCT, 1, sqrt(2), 1);
en = @(X) dct_energies(Y(X), N_DCT);

%N_try = 1;

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
        m = m + 1;
    end
end

%% Run on all files
M = length(f_name);
%inds = randperm(M);
%ssimVals = -1*ones(4, M);
%energies = zeros(length(f_name), N_DCT, N_DCT);
tic
t = 0;
prev_perc = 0;
cur_perc = 0;
cur_m = 42;
for m = cur_m:M
    try
        X_l = image_reader(f_name{m});
        evalc('[ssimVals(1, m),~,~, ~, ~] = dwt_opt_enc(X_l, 7, -1, 0.5);');
        evalc('[ssimVals(2, m),~,~,~] = jpegencdeclbtnlev(X_l, 4, 16, 0.5, sqrt(2), true, 16, 0.1, 0);'); 
        evalc('[ssimVals(3, m),~,~,~] = jpegencdeclbtnlev(X_l, 4, 16, 0.5, sqrt(2), true, 16, 0.2, 0);'); 
        evalc('[ssimVals(4, m),~,~,~] = jpegencdeclbtnlev(X_l, 4, 16, 0.5, sqrt(2), true, 16, 0.3, 0);'); 
    catch
        fprintf("Image %i failed. Skipping.\n", m);
        ssimVals(:,m) = -Inf;
    end
    %updates
    t_elapsed = toc;
    fprintf("Image %i\n", m);
    fprintf("Average time per compression: %0.2f minutes\n", t_elapsed/(60*(m-cur_m+1)));
    fprintf("Expected time left: %0.2f hours\n\n", (M - m)*t_elapsed/(60*60*(m-cur_m+1)));
    
    cur_perc = floor(100*m/M);
    if(cur_perc ~= prev_perc)
        prev_perc = cur_perc;
        fprintf("\n\n %i%% COMPLETE \n\n", prev_perc);
    end
end
%{
%% Get files from output
N_STORE = 1;
files = cell(N_DCT, N_DCT, 3, N_STORE);
ssimVals = zeros(2, N_DCT, N_DCT, 3, N_STORE);
for r = 1:N_DCT
    for c = 1:N_DCT
        [A, I] = sort(energies(:, r, c));
        m_ind = floor(0.5*(length(I) - N_STORE));
        
        for n = 1:N_STORE
            %highest energy
            files{r, c, 1, n} = f_name{I(1 + n)};
            X_l = image_reader(files{r, c, 1, n});
            [ssimVals(1, r, c, 1, n),~,~, ~, ~] = dwt_opt_enc(X_l, 7, -1, 1);
            evalc('[ssimVals(2, r, c, 1, n),~,~,~] = jpegencdeclbtnlev(X_l, 4, 16, 1, sqrt(2), true, 16, 0.2, 0);');
            
            
            %medium energy
            files{r, c, 2, n} = f_name{I(m_ind + n)};
            X_l = image_reader(files{r, c, 2, n});
            [ssimVals(1, r, c, 2, n),~,~, ~, ~] = dwt_opt_enc(X_l, 7, -1, 1);
            evalc('[ssimVals(2, r, c, 2, n),~,~,~] = jpegencdeclbtnlev(X_l, 4, 16, 1, sqrt(2), true, 16, 0.2, 0);');
            
            %high energy
            files{r,c,3,n} = f_name{I(length(I) - N_STORE + n)};
            X_l = image_reader(files{r,c,3,n});
            [ssimVals(1, r, c, 3, n),~,~, ~, ~] = dwt_opt_enc(X_l, 7, -1, 1);
            evalc('[ssimVals(2, r, c, 3, n),~,~,~] = jpegencdeclbtnlev(X_l, 4, 16, 1, sqrt(2), true, 16, 0.2, 0);');
        end
    end
end

save('test_results','ssimVals','energies','files')
%}    