function B = repmat_my(A,p1,p2,p3)

n = size(A);
if(length(n)==2 & p1==1 & p2==1 & p3>1)
    B = zeros(n(1),n(2),p3);
    for i = 1:p3
       B(:,:,i) = A; 
    end
end

if(length(n)==3 & p1==1 & p2>1 & p3==1)
    B = zeros(n(1),p2,n(3));
    for i = 1:p2
        B(:,i,:) = A;
    end
end

if(length(n)==3 & p1>1 & p2==1 & p3==1)
    B = zeros(p1,n(2),n(3));
    for i = 1:p1
        B(i,:,:) = A;
    end
end