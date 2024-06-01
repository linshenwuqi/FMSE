function results = FMSE_v(FCE_all, clsArr, clsNums, para_theta, Para_alpha)

if nargin < 5
    para_theta = 0.84;
    Para_alpha = 0.95;
end

% initial
M = length(clsArr);

FRI = computeFRI(FCE_all,clsArr, para_theta);

EFCA = computeEFCA(FCE_all, FRI, M, Para_alpha);

results = SCPP_v(EFCA, clsNums);
end