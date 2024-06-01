function [purity] = compute_Purity(labels_true, labels_pred)
clusters = unique(labels_pred);
labels_true = labels_true';
labels_pred = labels_pred';
labels_true = labels_true(:);
labels_pred = labels_pred(:);
count = [];

for c = 1:length(clusters)
	idx = find(labels_pred == c);
	temp = labels_true(idx);
	labels_tmp = reshape(temp,1,length(temp(:)));
	T=tabulate(labels_tmp);
	count = [count, max(T(:,2))];
end
purity = sum(count)/size(labels_true,1);
end
