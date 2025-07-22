function [ q_matrix_Alg, Lag_SR_mid] = Alg_WSR_PT_PowerAllocation( Weight_Total,QoS,eta,beta,theta,xi,rho,Psi,P_max ) 
global M
global N
global Channel

Thr_total = -1*ones(N,M);
Index_total = -1*ones(N,M);
Monot_total = -2*ones(N,M);
weight_SC = -1*ones(N,M);
alpha_SC = -1*ones(N,M);
xi_SC = -1*ones(N,M);
rho_SC = -1*ones(N,M);
Dev_Total = -1*ones(N+1,M);
Thr_Lambda = -1*ones(N+1,M);

for m=1:M
    [ Thr_SC, Index_begin, Monot_SC ] = Get_Thr_SC( Weight_Total(:,m), Channel(:,m), QoS(:,m), eta, xi, beta, theta, rho, Psi);
    Thr_total(1:length(Index_begin),m) = Thr_SC;
    Index_total(1:length(Index_begin),m) = Index_begin;
    Monot_total(1:length(Index_begin),m) = Monot_SC;
    for i = 1:length(Index_begin)
        n_begin = Index_begin(i);
        if i<length(Index_begin)
            n_end = Index_begin(i+1)-1;
        else
            n_end = N;
        end
        if Thr_total(i,m) > Psi(i,m)
            weight_SC(i,m) = Weight_Total(n_end,m);
            alpha_SC(i,m) = Channel(n_end,m);
            xi_SC(i,m) = xi(1,n_end,m);
            rho_SC(i,m) = rho(1,n_end,m);
        end
        if i>1
            Thr_Lambda(i,m) = (Thr_total(i,m)+rho(1,n_begin,m))/xi(1,n_begin,m);
        end
    end

    for i = 1:length(Index_begin)
        if i<length(Index_begin) && weight_SC(i,m)>0
            Dev_Total(i+1,m) = weight_SC(i,m)*alpha_SC(i,m)*xi_SC(i,m)/...
                (log(2)*(...
                1+alpha_SC(i,m)*xi_SC(i,m)*Thr_Lambda(i+1,m)-alpha_SC(i,m)*rho_SC(i,m)...
                ));
        end
        if i==length(Index_begin) && weight_SC(i,m)>0
            Dev_Total(i+1,m) = weight_SC(i,m)*alpha_SC(i,m)*xi_SC(i,m)/...
                (log(2)*(...
                1+alpha_SC(i,m)*xi_SC(i,m)*Psi(1,m)-alpha_SC(i,m)*rho_SC(i,m)...
                ));
        end
    end
end

Lag_SR_up = max(max(Weight_Total.*Channel));
Lag_SR_down = 0;
Lag_gap = 10^-8;

Index_t = zeros(M,1);
Power_Lag_save = zeros(M,1);

while Lag_SR_up-Lag_SR_down>Lag_gap
    Lag_SR_mid =(Lag_SR_up+Lag_SR_down)/2;
    Power_Lag = 0;
    for m=1:M
        Dev_count = 1;
        while Dev_count<N && Dev_Total(Dev_count+1,m)<Lag_SR_mid && Dev_Total(Dev_count+1,m)>0
            Dev_count = Dev_count+1;
        end
        weight_Lag = weight_SC(Dev_count,m);
        alpha_Lag = alpha_SC(Dev_count,m);
        xi_Lag = xi_SC(Dev_count,m);
        rho_Lag = rho_SC(Dev_count,m);
        Power_Lag_tmp = weight_Lag/(Lag_SR_mid*log(2))+rho_Lag/xi_Lag-1/(alpha_Lag*xi_Lag);
        Power_Lag_tmp = max(Power_Lag_tmp,Psi(1,m));
        Power_Lag = Power_Lag+Power_Lag_tmp;
        Power_Lag_save(m) = Power_Lag_tmp;
        Index_t(m) = Dev_count+1;
    end
    if Power_Lag>P_max
        Lag_SR_down = Lag_SR_mid;
    else
        Lag_SR_up = Lag_SR_mid;
    end

end

q_matrix_Alg = -1*ones(N,M);
q_matrix_Alg(1,:) = Power_Lag_save';

for m=1:M
    Index_Begin = Index_total(:,m);
    Monot_SC = Monot_total(:,m);
    for t_tmp = Index_t(m):sum(Index_Begin>0)
        n_begin = Index_Begin(t_tmp);
        if xi(1,n_begin,m)*q_matrix_Alg(1,m)-rho(1,n_begin,m)>=Thr_total(t_tmp,m)
            q_matrix_Alg(n_begin,m) = Thr_total(t_tmp,m);
        end
    end
    for n=1:N-1
        if q_matrix_Alg(n+1,m)<0
            q_matrix_Alg(n+1,m) = eta(n,m)*q_matrix_Alg(n,m)-beta(n,m);
        end
    end
end