function [ Thr_SC, Group_begin, Monot_list] = Get_Thr_SC( weight, alpha, QoS, eta_omit, xi_omit, beta_omit, theta_omit, rho_omit, Psi_omit)
global N;
eta = 2.^(-QoS);
beta = (1-eta)./alpha;
theta = (2.^QoS-1)./alpha;

xi = eye(N);
for n=N:-1:2
    for l=n:N
        xi(n-1,l) = eta(n-1)*xi(n,l);
    end
end

rho = zeros(N,N);
for n=N:-1:2
    for l=n:N
        rho(n-1,l) = rho(n,l) + xi(n,l)*beta(n-1);
    end
end

Psi = zeros( N,1 );
for n=1:N
    Psi(n) = (theta(N) + rho(n,N))/xi(n,N);
end

Thr_SC = zeros( N,1 );
Thr_SC(1) = Inf;
Monot_list = zeros( N,1 );
Monot_list(1) = 1;
Group_begin = 1:N;


for n=2:N
    [Monot_indicator, Thr_sub] = Group_Monot(weight, alpha, xi, rho, Psi, n, n);
    Thr_SC(n) = Thr_sub;
    Monot_list(n) = Monot_indicator;
end

Conv_indicator=0;
while (Conv_indicator==0)
    [Conv_indicator, Merge_index] = Conv_check(weight, alpha, xi, rho, Psi, Monot_list, Thr_SC, Group_begin);
    if Conv_indicator==0
        Index_begin = Group_begin(Merge_index-1);
        if Merge_index == length(Group_begin)
            Index_end = N;
        else
            Index_end = Group_begin(Merge_index+1)-1;
        end
        Thr_SC(Merge_index)=[];
        Monot_list(Merge_index)=[];
        Group_begin(Merge_index)=[];
        [Monot_indicator, Thr_sub] = Group_Monot(weight, alpha, xi, rho, Psi, Index_begin, Index_end);
        Thr_SC(Merge_index-1) = Thr_sub;
        Monot_list(Merge_index-1) = Monot_indicator;
    end
end

