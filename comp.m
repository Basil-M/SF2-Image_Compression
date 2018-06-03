%% FIRST FILE
clear all;
fname = 'lighthouse';

load(fname);
fprintf("Compressing %s to 40960 bits.\n", fname);
fname = fname(1:5);
comp_enc;
comp_dec;

Z_draw = Z;
%% SECOND FILE
clearvars -except Z_draw
fname = 'bridge';

load(fname);
fprintf("Compressing %s to 40960 bits.\n", fname);
fname = fname(1:5);
comp_enc;
comp_dec;

Z_draw = beside(Z_draw, Z);

%% THIRD FILE
clearvars -except Z_draw
fname = 'flamingo';

load(fname);
fprintf("Compressing %s to 40960 bits.\n", fname);
fname = fname(1:5);
comp_enc;
comp_dec;

Z_draw = beside(Z_draw, Z);

%% DRAW ALL
draw(Z_draw);