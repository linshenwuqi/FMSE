function [accuracy, ConMtx] = compute_ACC(result,gt)

n = length( gt );
clsNum = max(gt);
E = zeros( clsNum, clsNum );
for idx = 1 : n
    i1 = result(idx);
    i2 = gt(idx);
    E( i1, i2 ) = E( i1, i2 ) + 1;
end
ConMtx=E';
E=-E;
[C,~]=hungarian(E);
nMatch=0;
for i=1:clsNum
    nMatch=nMatch-E(C(i),i);
end

accuracy = nMatch/n;
