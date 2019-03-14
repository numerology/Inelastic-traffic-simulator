% Script to test the performance under a static setting
%% Params setup
nUsers = 6;
nSlices = 2;
nBasestations = 2;
capPerUser = ones(1, 6);
bsPerUser = [1 1 1 2 2 2];
shareDist = [1/2 1/2;1/2 1/2];
opBelongs = [1 1 2 1 2 2];
shareVec = [1 1];

%%
netSettings = [];
netSettings.users = nUsers;
netSettings.bsNS = nBasestations;
opSettings = [];
opSettings.s_o = shareVec;
opSettings.w_i = zeros(size(capPerUser));
opSettings.ops_belongs = opBelongs;
opSettings.shareDist = shareDist;

% nUsers = NetSettings.users;
% nSlices = size(OpSettings.s_o, 2);
% nBasestations = NetSettings.bsNS;
% shareVec = OpSettings.s_o;
% opBelongs = OpSettings.ops_belongs;
% remainingPerBs = ones(1, nBasestations);
% bidPerUser = zeros(size(OpSettings.w_i));
% shareDist = OpSettings.shareDist;

%%
[r,f,b] = flexibleSS(netSettings, opSettings, capPerUser, bsPerUser);
disp('Under static slicing')
f
b
[r,f,b] = biddingSCG(netSettings, opSettings, capPerUser, bsPerUser);
disp('Under share with guarantees')
f
b