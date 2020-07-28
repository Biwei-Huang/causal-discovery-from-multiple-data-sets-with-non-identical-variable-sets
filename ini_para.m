function initA = ini_para(Data,IDX)

M = length(Data);
n = size(IDX,1);
initA = zeros(n,n);
for m = 1:M
    X = Data{m};
    id = find(IDX(:,m)==1);
    
    c{m} = zeros(n,n);
    for i = 1:length(id)
        j = id(i);
        k = setdiff(id,j);
        for r = 1:length(k)
            p = k(r);
            a1 = lasso(X(p,:)',X(j,:)');
            c{m}(j,k) = a1(:,50)';
        end
    end
%     c{m} = (c{m}+c{m}')/2;
end

c_min = 999*(ones(n,n) - eye(n));
for i = 1:n
    for j = 1:n
        if(i~=j)
            for m = 1:M
                if(abs(c_min(i,j))>abs(c{m}(i,j)) & c{m}(i,j)~=0)
                    c_min(i,j) = c{m}(i,j);
                end
            end
            if(c_min(i,j)==999)
                c_min(i,j) = 0.1;
            end
        end
    end
end
initA = inv(eye(n)-c_min);
