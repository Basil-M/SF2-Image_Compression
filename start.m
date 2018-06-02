%% Simply adds all the paths to the repository and starts up clean with lighthouse

addpath('DWT');
addpath('COMMON');
addpath('IMAGE_TESTING');
addpath('LBT');
addpath('IMAGE_TESTING/IMAGES');

clear all;
load flamingo;
X_f = X;
load bridge;
X_b = X;
load lighthouse;
X_lh = X;
set(0,'defaultTextInterpreter','latex');
