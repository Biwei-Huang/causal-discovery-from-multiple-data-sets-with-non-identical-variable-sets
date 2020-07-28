% causal discovery from multiple datasets with non-identical variable sets
%     by maximizing the log-likelihood
% 4 variables, 3 data sets, 1000 sample size of each data set

clear all,clc,close all;
addpath(genpath(pwd))

%% generate data
rng(10);
M = 3; % total number of data sets
T = 1000*ones(1,M); % number of samples in each data set
n = 4; % total number of integrated variables
n0 = 3; % number of variables in each dataset
combs = nchoosek(1:n,n0);
tmp_idx = zeros(n,size(combs,1));
for i = 1:size(combs,1)
    tmp_idx(combs(i,:),i) = 1;
end

while(1)
    tmp = randperm(size(tmp_idx,2));
    IDX_tmp = tmp_idx(:,tmp(1:M));
    sign = zeros(1,M);
    for i = 1:M
        for j = i+1:M
            if(~isempty(intersect(IDX_tmp(:,i),IDX_tmp(:,j))))
                sign(i)=1;
                sign(j)=1;
            end
        end
    end
    if(sum(sign)==M & isempty(find(all(IDX_tmp==0,2))))
        IDX = IDX_tmp;
        break;
    end
end

mu = repmat([-0.5 0.5],[n,1]);
w = repmat([0.5 0.5],[n 1]);
sigma = repmat([0.1 0.1],[n 1]);

G0 = generate_structure(n);
B0 = (0.3*rand(n,n)+0.5).*G0';
A0 = inv(eye(n)-B0);
for m = 1:M
    E{m} = zeros(n, T(m));
    for i = 1:n
        MU = mu(i,:)';
        SIGMA = reshape(sigma(i,:),1,1,size(w,2));
        P = w(i,:);
        obj = gmdistribution(MU,SIGMA,P);
        
        E{m}(i,:) = random(obj, T(m));
    end
    I{m} = diag(IDX(:,m));
    x{m} = I{m}*A0*E{m};
end
Data = x;


%% causal discovery from multiple data sets with non-identifical sets of variables
% set the parameters
parsEM.thres = 1e-7;
parsEM.maxIter = 30000;
initNoise = 0.6;
parsEM.A = A0+initNoise*rand(n,n)-initNoise/2;
parsEM.mu = mu;
parsEM.sigma = sigma;
parsEM.w = w;
parsEM.updatePrior = 1;
parsEM.zeroMean = 0;
parsEM.noise = 1e-4;
parsEM.fast = 0;
tic;
filename = 'example1_result';
[A_hat, muHat, sigmaHat, wHat, loglAll] = ifaEM_my(Data, IDX, 2, parsEM, filename);
toc;


%% permutation to derive W and B matrices
W_hat = inv(A_hat);
% Try to permute the rows of W so that sum(1./abs(diag(W))) is minimized
fprintf('Performing row permutation...\n');
dims = n;
if dims<=8,
    fprintf('(Small dimensionality, using brute-force method.)\n');
    [Wp,rowp] = nzdiagbruteforce(W_hat);
else
    fprintf('(Using the Hungarian algorithm.)\n');
    [Wp,rowp] = nzdiaghungarian(W_hat);
end
fprintf('Done!\n');

% Divide each row of Wp by the diagonal element
estdisturbancestd = 1./diag(abs(Wp));
Wp = Wp./(diag(Wp)*ones(1,dims));

% Compute corresponding B
Best = eye(dims)-Wp; % estimated causal adjacency matrix over the integrated set of variables

save(filename,'Data','G0','B0','A0','IDX','Best');

