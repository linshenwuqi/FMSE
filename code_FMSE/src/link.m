function [sim, result,pdf] = link(sim, initial, init_k)

n = size(sim,1);
result = zeros(n,1);
pdf = zeros(n,init_k);
sqrt_n = ceil(sqrt(n));

for i = 1:init_k 
    result(initial(i)) = i;
    pdf(initial(i),i) = 1;
end

%% Processing unassigned samples: 1. Similarity ranking; 2. Merge into prototype sample set R;
while sum(result == 0) > 0
% -------------------------Initialize similarity matrix-----------------------------
    unlabel_idx = find(result == 0);
    labeled_num = n-length(unlabel_idx);
    sim_init = zeros(length(unlabel_idx),length(labeled_num));
    for i = 1:init_k
        cat_mark_i = result == i;
        unlabel_simi = sim(unlabel_idx,cat_mark_i); 
        sim_init(:,i) = max(unlabel_simi,[],2);
    end

    norm_sim_init = sim_init ./ sum(sim_init,2);
    re_sim_in = sort(sim_init,2,'descend');
    cha_sim_in = abs(re_sim_in(:,1) - re_sim_in(:,2));
    
    if sum(sum(cha_sim_in)) == 0
        new_sorted_idx = 1:1:length(cha_sim_in);
        result(unlabel_idx(new_sorted_idx)) = init_k+1;
    else
        [~,sorted_idx] = sort(cha_sim_in,'descend');
        new_sorted_idx = sorted_idx(1 : min(sqrt_n, length(sorted_idx)));
% ----------------------------Assigning samples--------------------------------------
        if isempty(new_sorted_idx)
            new_sorted_idx = 1:1:length(cha_sim_in);
        end
        [~,label] = max(sim_init(new_sorted_idx,:), [], 2);
        result(unlabel_idx(new_sorted_idx)) = label;
        pdf(unlabel_idx(new_sorted_idx),:) = norm_sim_init(new_sorted_idx,:);
    end
end