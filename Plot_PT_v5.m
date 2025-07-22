close all;
clear all;
clc;
tic;

load('PT_20250610.mat');

Power_dBm = Power_dB+30;
PT_dBm = Power_PT_dB+30;

PowerSum_dBm = PowerSum_record_PT;
PowerSum_dBm = 30+10*log10(PowerSum_dBm);
PT_Thr_dBm = 30+10*log10(PT_Thr_record);
PowerSum_Thr_dBm = 30+10*log10(PowerSum_record_PT_Thr);


PT_Thr_dBm_total = zeros(3,length(Power_dBm));
for p_dB = 1:length(Power_dBm)
    i=-P_T_gap_dB:P_T_gap_dB:P_T_gap_dB;
    i=i';
    PT_Thr_dBm_total(:,p_dB) = i+PT_Thr_dBm(p_dB);
end


PT_dBm_total = zeros(length(PT_dBm)+3,length(Power_dBm));
PowerSum_dBm_total = zeros(length(PT_dBm)+3,length(Power_dBm));
for p_dB = 1:length(Power_dBm)
    i_PT = 1;
    i_Thr = 1;
    i = 0;
    while i<length(PT_dBm)+3
        i = i+1;
        if i_Thr>3 || PT_dBm(i_PT)<PT_Thr_dBm_total(i_Thr,p_dB)            
            PT_dBm_total(i,p_dB) = PT_dBm(i_PT);
            PowerSum_dBm_total(i,p_dB) = PowerSum_dBm(i_PT,p_dB);
            i_PT = i_PT+1;
        else            
            PT_dBm_total(i,p_dB) = PT_Thr_dBm_total(i_Thr,p_dB);
            PowerSum_dBm_total(i,p_dB) = PowerSum_Thr_dBm(i_Thr,p_dB);
            i_Thr = i_Thr+1;
        end
    end
end

figure1=figure;
plot(PT_dBm_total(:,1),PowerSum_dBm_total(:,1),':p','Color',[149/255 44/255 149/255],'LineWidth',2,'MarkerSize',6);
hold on;
plot(PT_dBm_total(:,2),PowerSum_dBm_total(:,2),':p','Color',[51/255 161/255 201/255],'LineWidth',2,'MarkerSize',6);
plot(PT_dBm_total(:,3),PowerSum_dBm_total(:,3),':p','Color',[149/255 149/255 44/255],'LineWidth',2,'MarkerSize',6);

grid on;
box on;
plot([PT_Thr_dBm_total(2,1),PT_Thr_dBm_total(2,1)],[PowerSum_Thr_dBm(2,1)-1.5,PowerSum_Thr_dBm(2,1)+1.5],'-.','Color',[149/255 44/255 149/255],'LineWidth',2);
plot([PT_Thr_dBm_total(2,2),PT_Thr_dBm_total(2,2)],[PowerSum_Thr_dBm(2,2)-1.5,PowerSum_Thr_dBm(2,2)+1.5],'-.','Color',[51/255 161/255 201/255],'LineWidth',2);
plot([PT_Thr_dBm_total(2,3),PT_Thr_dBm_total(2,3)],[PowerSum_Thr_dBm(2,3)-1.5,PowerSum_Thr_dBm(2,3)+1.5],'-.','Color',[149/255 149/255 44/255],'LineWidth',2);


legend1=legend(["$\mathit{P_{max}}$=30dBm",...
    "$\mathit{P_{max}}$=35dBm",...
    "$\mathit{P_{max}}$=40dBm" ...
    "$\mathit{P_{T}}$ Threshold",...
    "$\mathit{P_{T}}$ Threshold",...
    "$\mathit{P_{T}}$ Threshold"]);
legend1.NumColumns = 2;
legend1.ItemTokenSize=[16,12];

set(legend1,'Position',[0.35 0.15 0.43 0.15],'Interpreter','latex');
hold off;
xlabel('Circuit Power (dBm)');
ylabel('Total Transmit Power (dBm)');
axis([30,55,24,42]);
FIGNAME = 'main_PT_20250610';
PrintFigToPaper('-depsc', FIGNAME, 16, 'Times New Roman', 7, 1, 0);