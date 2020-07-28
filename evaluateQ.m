function [posq, marginalx] = evaluateQ(Y, A, mu, sigma, Dim, dimG, numGauss, qidMat, priorq, Lambda)

N = size(Y,2);
Dim1 = Dim - dimG;
condProb = zeros(N, numGauss^Dim1);
for q = 1:size(condProb,2)
    indicies = sub2ind(size(mu),1:Dim,qidMat(q,:));
    muq = mu(indicies);
    sigmaq = sigma(indicies);
    muq = repmat(A * muq',[1 N]);
    U = A * diag(sigmaq.^2) * A' + Lambda;
    condProb(:,q) = (2*pi)^(-size(Y,1)*0.5) * det(U)^-0.5 * exp(-0.5*diag((Y-muq)'*inv(U)*(Y-muq)));
end
marginalx = sum(condProb.*repmat(priorq,[N,1]),2);
posq = condProb .* repmat(priorq,[N,1]) ./ repmat((marginalx+eps^20),[1 size(qidMat,1)]);