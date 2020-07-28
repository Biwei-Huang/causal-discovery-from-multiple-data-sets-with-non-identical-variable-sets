function priorq = priorLatentQ(w, Dim, dimG, numGauss, qidMat)
% Calculate the prior distribution of q
Dim1 = Dim - dimG;
priorq = zeros(1, numGauss^Dim1);
for q = 1:length(priorq)
    indicies = sub2ind(size(w),1:Dim,qidMat(q,:));
    wq = w(indicies);
    priorq(q) = prod(wq);
end