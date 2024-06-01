function [lden] = local_dens(dis, dens)

n = length(dens);
lden = zeros(n,1);
maxdis = max(max(dis)) + 1;

index_maxdens = find(dens == max(dens));
lden(index_maxdens(1)) = maxdis + 1;

if length(index_maxdens)>1
    dens(index_maxdens(2:end)) = dens(index_maxdens(2:end)) - 0.0001;
end

for i = [1:index_maxdens(1)-1, index_maxdens(1)+1:n]
    now_dis = dis(i,:);
    now_dis(i) = maxdis + 1;
    lo_deni = dens == dens(i);
    
    if sum(lo_deni) > 1
        dens(lo_deni) = dens(lo_deni)-0.000001;
        dens(i) = dens(i)+0.000001;
    end
    
    index = dens > dens(i);
    no_dis = min(now_dis(index));
    lden(i) = no_dis(1);
end