function [ Rate_Total ] = Rate_Total( q_matrix )       
global M
global N
global Channel
global Bandwidth

Rate_Total = zeros(N,M);
for m=1:M    
    for n=1:N-1
        Rate_Total(n,m) = log2(1+Channel(n,m)*q_matrix(n,m))-log2(1+Channel(n,m)*q_matrix(n+1,m));
    end
    Rate_Total(N,m) = log2(1+Channel(N,m)*q_matrix(N,m));
end
Rate_Total = Bandwidth*Rate_Total;

