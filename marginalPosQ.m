function marginalPosq = marginalPosQ(posq, Dim, dimG, numGauss, qidMat)
% Marginal posterior of q
T = size(posq, 1);
marginalPosq = zeros(Dim, numGauss, T);
for t = 1:T
    for i = 1:size(marginalPosq,1)
        for j = 1:size(marginalPosq,2)
            posqt = posq(t,:);
            ind = qidMat(:,i)==j;
            marginalPosq(i,j,t) = sum(posqt(ind));
        end
    end
end
end