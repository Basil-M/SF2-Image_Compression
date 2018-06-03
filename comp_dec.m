
%load(fname+'cmp');
load(strcat(fname,'cmp'));
%print number of bits
n_bits = sum(vlc(:,2));
fprintf('Number of bits for vlc:\t\t\t\t\t %i\n', n_bits);
if ~exist('rise','var')
    N = 4;
    M = 16;
    rise1 = 1;
    s = sqrt(2);
    dcbits = 16;
    ratio2 = 0;
    %LBT DECODE
    evalc('Z = jpegdeclbtnlev(vlc, q_opt_lbt, N, M, rise1, s, bits, huffval, dcbits, ratio,ratio2);')
    
    % Header information
        % ratio - float - 32 bits
        % q_opt - float - 32 bits
        % bits/huffval  - 1424 bits
    n_bits = n_bits + 32 + 32 + 1424;
    fprintf('Number of bits for header information:\t %i\n',32+32 + 1424);
else
    %% Decoder has variables
        % q_opt
        % rise
        % vlc
        % bits
        % huffval
    evalc('Z = dwt_dec(vlc, 7, -1, q_opt, rise, 0, bits, huffval);');
    
    % Header information
        % q_opt - float - 32 bits
        % rise  - float - 32 bits
        % bits/huffval  - 1424 bits
    n_bits = n_bits + 32 + 32 + 1424;
    fprintf('Number of bits for header information:\t %i\n',32+32 + 1424);
end
% print total number of bits
fprintf('Total number of bits: \t\t\t\t\t %i\n',n_bits);
%print RMS error (same for both)
fprintf('RMS error: \t\t\t %0.2f\n',std(Z(:) - X(:)));
