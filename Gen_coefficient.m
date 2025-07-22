function [ eta,beta,theta,xi,rho,Psi ] = Gen_coefficient( QoS )        
global M
global N
global Channel

eta = 2.^(-QoS);
beta = (1-eta)./Channel;
theta = (2.^QoS-1)./Channel;
xi = zeros(N,N,M);
rho = zeros(N,N,M);
Psi = zeros(N,M);


for m=1:M
    xi(:,:,m) = eye(N);
    for n=N:-1:2
        for l=n:N
            xi(n-1,l,m) = eta(n-1,m)*xi(n,l,m);
        end
    end
    rho(:,:,m) = zeros(N,N);
    for n=N:-1:2
        for l=n:N
            rho(n-1,l,m) = rho(n,l,m) + xi(n,l,m)*beta(n-1,m);
        end
    end
    for n=1:N
        Psi(n,m) = (theta(N,m) + rho(n,N,m))/xi(n,N,m);
    end
end