function [ q_matrix ] = Q_Alloc( Power_Lag_save,Index_total,Monot_total,Index_t,Thr_total,eta,beta,xi,rho )        %第4条曲线的目标函数,此处q不包含q1
global M
global N

q_matrix = -1*ones(N,M);
q_matrix(1,:) = Power_Lag_save;


for m=1:M
    Index_Begin = Index_total(:,m);
    Monot_SC = Monot_total(:,m);
    for t_tmp = Index_t(m):sum(Index_Begin>0)
        n_begin = Index_Begin(t_tmp);
        if xi(1,n_begin,m)*q_matrix(1,m)-rho(1,n_begin,m)>=Thr_total(t_tmp,m)
            q_matrix(n_begin,m) = Thr_total(t_tmp,m);
        end
    end
    for n=1:N-1
        if q_matrix(n+1,m)<0
            q_matrix(n+1,m) = eta(n,m)*q_matrix(n,m)-beta(n,m);
        end
    end
end