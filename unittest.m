%% Test of bestresponseiteration
clc, clear all
% case 1:
cBid = 1/4 * ones(1, 4);
opBelongs = [1 1 2 2];
bs = [1 2 1 2];
shareDist = 1 / 2 * ones(2, 2);
shareVec = [1 1];
capacityPerUser = ones(1, 4);

nextBid1 = bestresponseiteration(cBid, 1, shareDist, shareVec, opBelongs, bs,...
    capacityPerUser);
nextBid2 = bestresponseiteration(nextBid1, 2, shareDist, shareVec, opBelongs, ...
    bs, capacityPerUser);
expectedBid1 = [0.5 0.5 0.25 0.25];
expectedBid2 = [0.5 0.5 0.5 0.5];
assert(norm(nextBid1 - expectedBid1) < 1e-5 && ...
    norm(nextBid2 - expectedBid2) < 1e-5, ...
    'unexpected output in test case 1');




disp('All case finished!');