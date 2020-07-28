function [posMeanE, posCovE, posMeanEqi, posCovEqi, posCovEq0] = evaluateE_my(Y, posq, A, Lambda, ...
    mu, sigma, Dim, dimG, numGauss, qidMat)

N = size(Y,2);
Dim1 = Dim - dimG;
posMeanE = zeros(Dim, N, numGauss^Dim1);
posCovEq = zeros(Dim, Dim, N, numGauss^Dim1);
posCovEq0 = zeros(Dim, Dim, numGauss^Dim1);
posMeanEqi = zeros(Dim, numGauss, N);
posCovEqi = zeros(Dim, numGauss, N);

for q = 1:size(posMeanE,3)
    indicies = sub2ind(size(mu),1:Dim,qidMat(q,:));
    muq = mu(indicies);
    sigmaq = sigma(indicies);
    sigmaq = sigmaq.^2;
    muq = repmat(muq',[1, N]);
    posMeanE(:,:,q) = muq + diag(sigmaq) * A' * ((A*diag(sigmaq)*A'+Lambda) \ (Y-A*muq));
    posCovEq0(:,:,q) = (diag(sigmaq) - diag(sigmaq)*A'*...
        ((A*diag(sigmaq)*A'+Lambda) \ A * diag(sigmaq)));
    rshPosMeanEV = reshape(posMeanE(:,:,q),Dim,1,N);
    rshPosMeanEH = reshape(posMeanE(:,:,q),1,Dim,N);

    posCovEq(:,:,:,q) = repmat_my(posCovEq0(:,:,q),1,1,N)+repmat_my(rshPosMeanEV,1,Dim,1).*...
    repmat_my(rshPosMeanEH,Dim,1,1);
    posMeanE(:,:,q) = repmat(posq(:,q)',[Dim 1]) .* posMeanE(:,:,q);
end
for q = 1:numGauss
    indicies = qidMat==q;
    indicies = reshape(indicies', Dim, 1, numGauss^Dim1);
    indicies = repmat_my(indicies, 1, N, 1);
    posMeanEqi0 = sum(indicies.*posMeanE, 3);
    posMeanEqi(:,q,:) = reshape(posMeanEqi0, Dim, 1, N);
end
posMeanE = sum(posMeanE,3);

posCovEqDiag = posCovEq(logical(repmat(eye(size(posCovEq(:,:,1))),[1 1 size(posCovEq,3)...
size(posCovEq,4)])));
posCovEqDiag = reshape(posCovEqDiag, [Dim N numGauss^Dim1]);
posqRep = repmat(reshape(posq,[1 N numGauss^Dim1]),[Dim 1 1]);
for q = 1:numGauss
    indicies = qidMat==q;
    indicies = reshape(indicies', Dim, 1, numGauss^Dim1);
    indicies = repmat_my(indicies, 1, N, 1);
    posCovEqi0 = sum(indicies.*posqRep.*posCovEqDiag, 3);
    posCovEqi(:,q,:) = reshape(posCovEqi0, Dim, 1, N);
end
   
for q = 1:size(posCovEq,4)
    posCovEq(:,:,:,q) = repmat(reshape(posq(:,q),[1 1 N]), [Dim Dim 1]) .* posCovEq(:,:,:,q);
end
posCovE = sum(posCovEq,4);