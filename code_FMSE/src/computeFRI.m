function FRI = computeFRI(FCE_all, clsArr, para_theta)

if nargin < 3
    para_theta = 0.84;
end

% uncertainty of the sample with respect to the cluster
[N, ~] = size(FCE_all);
M = length(clsArr);

M_entro = zeros(N,M);
startCol = 1;

for i = 1:M
    numCols = clsArr(i);
    endCol = startCol + numCols - 1;
    
    FCE_i = FCE_all(:, startCol:endCol);
    M_entro(:,i) = -sum(FCE_i .* log2(FCE_i), 2);

    startCol = endCol + 1;
end

% M_entro(FCE_all)
ETs = sum(M_entro, 2);
FRI = exp(-ETs ./ para_theta ./ (M-1));

end