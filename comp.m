%% ADD PATHS
addpath('DWT');
addpath('COMMON');
addpath('LBT');

%% FILENAMES
fnames = cell(1,3);
fnames{1} = 'bridge';
fnames{2} = 'flamingo';
fnames{3} = 'competition';

for cur_f = 2:-1:1
    clearvars -except fnames cur_f;
    fname = fnames{cur_f};

    load(fname);
    X = double(X);
    fprintf('Compressing %s to 40960 bits.\n', fname);
    fname = fname(1:5);
    comp_enc;
    clearvars -except fnames fname cur_f;
    comp_dec;
end

%% LOAD IMAGES
%Upload your competition  entry as a single .mat file, 
%named  GroupN.mat -- where N is your group letter. 
%The .mat file should contain three matrices, 
%X (the decoded competition image), X1 (the decoded Bridge image)
%and X2 (the decoded Flamingo image). All images should remain 
% of size 256x256 with values   in the range 0 to 255.

clearvars -except fnames;
fname = fnames{1};
load(strcat(fname(1:5),'dec'));
X1 = Z;

fname = fnames{2};
load(strcat(fname(1:5),'dec'));
X2 = Z;

fname = fnames{3};
load(strcat(fname(1:5),'dec'));
X = Z;
clearvars -except X X1 X2;
save('Group9');