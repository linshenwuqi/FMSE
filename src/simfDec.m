function simf_approx = simfDec(FCA, Para_alpha)

if nargin < 2
    Para_alpha = 0.95;
end

[U, S, V] = mySVD(FCA);

singular_values = diag(S);
energy = cumsum(singular_values.^2) / sum(singular_values.^2);
k = find(energy >= Para_alpha, 1);

Uk = U(:, 1:k);
Sk = S(1:k, 1:k);
Vk = V(:, 1:k);

simf_approx = Uk * Sk * Vk';
end