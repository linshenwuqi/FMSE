function EFCA = computeEFCA(FCE_all, FRI, M, Para_alpha)

% inital Fuzzy Co-Association Martix
simf = (FCE_all * FCE_all') ./ M;

% SVD Reconstruct FCA Martix
simf_approx = simfDec(simf, Para_alpha);

% Enhanced Fuzzy Co-Association Martix
EFCA = ((FRI' .* simf_approx) + (FRI' .* simf_approx)') / 2;
end