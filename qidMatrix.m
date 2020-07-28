function qidMat = qidMatrix(numGauss, Dim, dimG)
% Get all the qs
Dim = Dim - dimG;
qidMat = zeros(numGauss^(Dim), Dim);
for i = size(qidMat,2):-1:1
    id = size(qidMat,2) - i;
    qidMatTemp = [];
    for j = 1:numGauss
        qidMatTemp = [qidMatTemp; j*ones(numGauss^id,1)];
    end
    qidMat(:,i) = repmat(qidMatTemp, [numGauss^(Dim-id-1),1]);
end
qidMat = [qidMat,ones(numGauss^Dim, dimG)];
end