function G = generate_structure(N)

p = 0.35;

G = triu(ones(N),1).*binornd(1,p,N,N);

for i = 2:N

    j = randi([1 i-1]);

    G(j,i) = 1;

end
    
end
