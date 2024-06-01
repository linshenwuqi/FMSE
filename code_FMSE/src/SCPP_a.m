function [results, res, init, pdf] = SCPP_a(simf, k)
% Inputs:
%   - simf: self-enhanced similarity matrix
%   - k: number of clusters
% Outputs:
%   - results: final clustering result
%   - res: cluster assignments after initial clustering
%   - init: Prototype sample set
%   - pdf: cluster representation matrix

n = size(simf, 1);

%% Density Peak
dens = diag(simf);
[lden] = local_dens(1 - simf, dens); % Calculate local density
use_lden = find(dens >= mean(dens));
[~, use_initial] = sort(lden(use_lden), 'descend');
init = use_lden(use_initial);

init_k = floor(sqrt(n));
if init_k < k
    init_k = k;
end
init = init(1:init_k); % Choose the top initial centroids

%% Assign samples to clusters and get initial clustering result
[simf, res, pdf1] = link(simf, init, init_k);

%% Single-linkage Hierarchical Clustering - Merge clusters
if init_k == k
    results = res;
else
    simij = zeros(1, (init_k - 1) * (init_k - 2) / 2 + (init_k - 1));
    for i = 1:init_k
        locat_i = res == i;
        for j = i + 1:init_k
            locat_j = res == j;
            sim_ij = simf(locat_i, locat_j);
            ij = sum(sum(sim_ij)) / (sum(locat_i) * sum(locat_j));
            simij((init_k - 1) * (i - 1) - ((i - 1) * (i - 2) / 2) + j - i) = ij;
        end
    end
    % Construct distance vector dis
    simij = simij ./ max(simij);
    dis = 1 - simij;
    % Construct tree structure and perform clustering
    Z = linkage(dis, 'single');
    clu = cluster(Z, k);
    % Assign each data point to the corresponding cluster center and calculate the pdf value for each cluster
    CC = zeros(n, 1);
    pdf = zeros(n, k);
    for i = 1:init_k
        biaoji = clu(i);
        M = res == i;
        CC(M, :) = biaoji;
    end
    for i = 1:k
        biaoji = clu == i;
        pdf(:, i) = sum(pdf1(:, biaoji), 2);
    end
    results = CC;
end
end
