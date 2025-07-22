clear;
clc;
tic;
%Created by Qian Peng in Southeast University in 2025 Jun. 10


Power_dB = 0:5:10;                           %power budget in dB

Power_PT_dB = 0:5:25;                %constant power in dB
Power_PT_dB = Power_PT_dB';
Level_PT = length(Power_PT_dB);
Level_Power = length(Power_dB);
Power_PT = zeros(Level_PT,1);
Power_budget = 10.^(Power_dB/10);

for lv_PT=1:Level_PT
    Power_PT(lv_PT) = 10^(Power_PT_dB(lv_PT)/10);
end

global N
global M
global P_T
global Bandwidth

N = 3;                                          %Number of UE per channel
M = 4;                                          %Number of channels
P_T = 1;
P_T_gap_dB = 0.5;
Radius_max = 500;                                   %Maximum Radius of cell
Radius_min = 20;                                    %Minimum Radius of cell
Bandwidth = 1;

Distance_record_PT = zeros(N,M);        
Channel_record_PT = zeros(N,M);      
Qower_record_PT = zeros(N,M,Level_PT,Level_Power);  
Rate_record_PT = zeros(N,M,Level_PT,Level_Power);   

PowerSum_record_PT = zeros(Level_PT,Level_Power);     
WEE_record_PT = zeros(Level_PT,Level_Power);        
Lag_SR_record = zeros(Level_Power,1);
PT_Thr_record = zeros(Level_Power,1);

Qower_record_PT_Thr = zeros(N,M,3,Level_Power);  
PowerSum_record_PT_Thr = zeros(3,Level_Power);     
WEE_record_PT_Thr = zeros(3,Level_Power);       


global Channel

[Channel,Distance]=Channel_coefficient(N,M,Radius_max,Radius_min,Bandwidth);
Channel_record_PT=Channel;

Weight_Rnd = Gen_Weight_Rnd( Channel );
QoS_Rnd = 2*rand(N,M);
[ eta_Rnd,beta_Rnd,theta_Rnd,xi_Rnd,rho_Rnd,Psi_Rnd ] = Gen_coefficient( QoS_Rnd );

for lv_PT=1:Level_PT
    for lv=1:Level_Power
        P_T = Power_PT(lv_PT);
        %Optimal Solution to WEE Problem
        [ q_matrix_WEE_Alg_Rnd,G_count_avg,D_count ] = Alg_WEE_PowerAllocation( Weight_Rnd,QoS_Rnd,eta_Rnd,beta_Rnd,theta_Rnd,xi_Rnd,rho_Rnd,Psi_Rnd,Power_budget(lv));
        Rate_WEE_Alg_Rnd = Rate_Stat( q_matrix_WEE_Alg_Rnd);
        Qower_record_PT(:,:,lv_PT,lv) = q_matrix_WEE_Alg_Rnd;
        PowerSum_record_PT(lv_PT,lv) = sum(q_matrix_WEE_Alg_Rnd(1,:));
        Rate_record_PT(:,:,lv_PT,lv) = Rate_WEE_Alg_Rnd;
        WEE_record_PT(lv_PT,lv) = sum(sum(Rate_WEE_Alg_Rnd.*Weight_Rnd))/(PowerSum_record_PT(lv_PT,1)+P_T);
    end
end
for lv=1:Level_Power
    [ q_matrix_WSR_Alg_Rnd, Lag_WSR_Rnd] = Alg_WSR_PT_PowerAllocation( Weight_Rnd,QoS_Rnd,eta_Rnd,beta_Rnd,theta_Rnd,xi_Rnd,rho_Rnd,Psi_Rnd,Power_budget(lv));
    Rate_WSR_Alg_Rnd = Rate_Stat( q_matrix_WSR_Alg_Rnd);
    Lag_SR_record(lv) = Lag_WSR_Rnd;
    PT_Thr_record(lv) = sum(sum(Rate_WSR_Alg_Rnd.*Weight_Rnd))/Lag_WSR_Rnd-Power_budget(lv);
end

for lv=1:Level_Power
    for lv_PT=1:3
        P_T = 10^(((lv_PT-2)*P_T_gap_dB)/10)*PT_Thr_record(lv);
        %Optimal Solution to WEE Problem
        [ q_matrix_WEE_Alg_Rnd,G_count_avg,D_count ] = Alg_WEE_PowerAllocation( Weight_Rnd,QoS_Rnd,eta_Rnd,beta_Rnd,theta_Rnd,xi_Rnd,rho_Rnd,Psi_Rnd,Power_budget(lv));
        Rate_WEE_Alg_Rnd = Rate_Stat( q_matrix_WEE_Alg_Rnd);
        Qower_record_PT_Thr(:,:,lv_PT,lv) = q_matrix_WEE_Alg_Rnd;
        PowerSum_record_PT_Thr(lv_PT,lv) = sum(q_matrix_WEE_Alg_Rnd(1,:));
        WEE_record_PT_Thr(lv_PT,lv) = sum(sum(Rate_WEE_Alg_Rnd.*Weight_Rnd))/(PowerSum_record_PT_Thr(lv_PT,1)+P_T);
    end
end
PT_Thr_record_dBm = 30+10*log10(PT_Thr_record);
PT_Thr_gap = zeros(3,1);
for lv=1:Level_Power
    PT_Thr_gap(lv) = PT_Thr_record_dBm(lv)-5*floor(PT_Thr_record_dBm(lv)/5);
end
toc;

if ~exist('PT_20250610.mat','file')
    save('PT_20250610');
    disp('file saved');
else
    disp('The file already exists, please confirm');
end