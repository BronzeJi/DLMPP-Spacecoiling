%% 
% partition top radius r3, 
% the neck diameter dt 
% the neck length L 
clear
close all

dB=@(a)(20*log10(a/(2*10^(-5))));

%environment constant
envir_cons.c=343;%speed of sound in m/s
envir_cons.mu=1.84*10^(-5);%mue pa*s
envir_cons.rho_c=1.213;% air density

%frequency
f=linspace(100,3000,5000);
deltaf=(3000-100)/5000;
%target frequency
f_target=linspace(10,2000,5000);



%%
figure("Name","Absorption_t");
[Plate1,Plate2]=set_default;
alpha=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha,"LineWidth",5);hold on;
area=deltaf*alpha;

Plate1.t=0.7*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha);

Plate1.t=0.4*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,":","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

Plate1.t=1.3*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-.","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

Plate1.t=1.6*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"--","LineWidth",5);hold off;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

ylim([0,1.1]);
xticks(200:300:3000);  % Specify desired tick marks
set(gca, 'XTickLabel', arrayfun(@(v) sprintf('%g', v/1000), xticks, 'UniformOutput', false));
xlabel("Frequency in kHz");ylabel("Absorption coefficient");legend("1mm","0.7mm","0.4mm","1.3mm","1.6mm","show","Orientation","horizontal","Location","southoutside");
box on
legend box off
fontsize(gcf,20,"points");
set(gcf,'position',[10,10,600,400])
saveas(gcf,'Absorption_t.png');


degree_t=0.3/1;
rate_Absorption_t=change_area./(4*degree_t*area);

figure(10);
plot(f,rate_Absorption_t,"LineWidth",4);
grid on
hold off

rate_Absorption_t_1=change_area/(4*degree_t);
figure(11);
plot(f,rate_Absorption_t_1,"LineWidth",4);
grid on
hold off
clear change_area

%%
figure("Name","Absorption_R3");
[Plate1,Plate2]=set_default;
alpha=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha,"LineWidth",5);hold on;
area=deltaf*alpha;

Plate2.R3=4*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha);


Plate2.R3=3*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,":","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;


Plate2.R3=6*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-.","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;


Plate2.R3=7*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"--","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

ylim([0,1.1]);
xticks(200:300:3000);  % Specify desired tick marks
set(gca, 'XTickLabel', arrayfun(@(v) sprintf('%g', v/1000), xticks, 'UniformOutput', false));
xlabel("Frequency in kHz");ylabel("Absorption coefficient");legend("5mm","4mm","3mm","6mm","7mm","show","Orientation","horizontal","Location","southoutside");
box on
legend box off
fontsize(gcf,20,"points");
set(gcf,'position',[10,10,600,400])
saveas(gcf,'Absorption_R3.png');


degree_R3=1/5;
rate_Absorption_R3=change_area./(4*degree_R3*area);
figure(10);
hold on
plot(f,rate_Absorption_R3,"LineWidth",4);
hold off

rate_Absorption_R3_1=change_area/(4*degree_R3);
figure(11);
hold on
plot(f,rate_Absorption_R3_1,"LineWidth",4);
hold off
%thickenall_big
clear change_area

%%
figure("Name","Absorption_d3");
[Plate1,Plate2]=set_default;
alpha=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha,"LineWidth",5);hold on;
area=deltaf*alpha;

Plate2.d3=0.3*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha);


Plate2.d3=0.1*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,":","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

Plate2.d3=0.7*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-.","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

Plate2.d3=0.9*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"--","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

ylim([0,1.1]);
xticks(200:300:3000);  % Specify desired tick marks
set(gca, 'XTickLabel', arrayfun(@(v) sprintf('%g', v/1000), xticks, 'UniformOutput', false));
xlabel("Frequency in kHz");ylabel("Absorption coefficient");legend("0.5mm","0.4mm","0.3mm","0.6mm","0.7mm","show","Orientation","horizontal","Location","southoutside");
box on
legend box off
fontsize(gcf,20,"points");
set(gcf,'position',[10,10,600,400])
saveas(gcf,'Absorption_d3.png');


degree_d3=0.2/0.5;
rate_Absorption_d3=change_area./(4*degree_d3*area);

figure(10);
hold on
plot(f,rate_Absorption_d3,"LineWidth",5);
hold off

rate_Absorption_d3_1=change_area/(4*degree_d3);
figure(11);
hold on
plot(f,rate_Absorption_d3_1,"LineWidth",5);
hold off
clear change_area
%thickenall_big


%%
figure("Name","Absorption_d2");
[Plate1,Plate2]=set_default;
alpha=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha,"LineWidth",5);hold on;

Plate2.d2=0.4*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha);

Plate2.d2=0.3*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,":","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

Plate2.d2=0.6*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-.","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

Plate2.d2=0.7*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"--","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

ylim([0,1.1]);
xticks(200:300:3000);  % Specify desired tick marks
set(gca, 'XTickLabel', arrayfun(@(v) sprintf('%g', v/1000), xticks, 'UniformOutput', false));
xlabel("Frequency in kHz");ylabel("Absorption coefficient");legend("0.5mm","0.4mm","0.3mm","0.6mm","0.7mm","show","Orientation","horizontal","Location","southoutside");
box on
legend box off
fontsize(gcf,20,"points");
set(gcf,'position',[10,10,600,400])

degree_d2=0.1/0.5;
rate_Absorption_d2=change_area./(4*degree_d2*area);
saveas(gcf,'Absorption_d2.png');

figure(10);
hold on
plot(f,rate_Absorption_d2,"-","LineWidth",5);
hold off

rate_Absorption_d2_1=change_area/(4*degree_d2);
figure(11);
hold on
plot(f,rate_Absorption_d2_1,"LineWidth",5);
hold off
%thickenall_big
clear change_area




%%
figure("Name","Absorption_p2");
[Plate1,Plate2]=set_default;
alpha=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha,"LineWidth",5);hold on;

Plate2.p2=5*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha);

Plate2.p2=15*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,":","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

Plate2.p2=20*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-.","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

Plate2.p2=25*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"--","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

ylim([0,1.1]);
xticks(200:300:3000);  % Specify desired tick marks
set(gca, 'XTickLabel', arrayfun(@(v) sprintf('%g', v/1000), xticks, 'UniformOutput', false));
xlabel("Frequency in kHz");ylabel("Absorption coefficient");legend("1%","0.5%","1.5%","2%","2.5%","show","Orientation","horizontal","Location","southoutside");
box on
legend box off
fontsize(gcf,20,"points");
set(gcf,'position',[10,10,600,400])
saveas(gcf,'Absorption_p2.png');


degree_p2=0.005/0.01;
rate_Absorption_p2=change_area./(4*degree_p2*area);


figure(10);
hold on
plot(f,rate_Absorption_p2,"--","LineWidth",5);
hold off

rate_Absorption_p2_1=change_area/(4*degree_p2);
figure(11);
hold on
plot(f,rate_Absorption_p2_1,"LineWidth",5);
hold off
%thickenall_big
clear change_area

%%
figure("Name","Absorption_p3");
[Plate1,Plate2]=set_default;
alpha=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha,"LineWidth",5);hold on;

Plate2.p3=5*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha);

Plate2.p3=15*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,":","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

Plate2.p3=20*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-.","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

Plate2.p3=25*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"--","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

ylim([0,1.1]);
xticks(200:300:3000);  % Specify desired tick marks
set(gca, 'XTickLabel', arrayfun(@(v) sprintf('%g', v/1000), xticks, 'UniformOutput', false));
xlabel("Frequency in kHz");ylabel("Absorption coefficient");legend("1%","0.5%","1.5%","2%","2.5%","show","Orientation","horizontal","Location","southoutside");
box on
legend box off
fontsize(gcf,20,"points");
set(gcf,'position',[10,10,600,400])
saveas(gcf,'Absorption_p3.png');



degree_p3=0.005/0.01;
rate_Absorption_p3=change_area./(4*degree_p3*area);


figure(10);
hold on
plot(f,rate_Absorption_p3,"-.","LineWidth",5);
hold off

rate_Absorption_p3_1=change_area/(4*degree_p3);
figure(11);
hold on
plot(f,rate_Absorption_p3_1,"LineWidth",5);
hold off
%thickenall_big
clear change_area

%%
figure("Name","Absorption_D2");
[Plate1,Plate2]=set_default;
alpha=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha,"LineWidth",5);hold on;

Plate2.D2=15*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha);

Plate2.D2=10*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,":","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

Plate2.D2=25*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-.","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

Plate2.D2=30*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"--","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;


ylim([0,1.1]);
xticks(200:300:3000);  % Specify desired tick marks
set(gca, 'XTickLabel', arrayfun(@(v) sprintf('%g', v/1000), xticks, 'UniformOutput', false));
xlabel("Frequency in kHz");ylabel("Absorption coefficient");legend("20mm","15mm","10mm","25mm","30mm","show","Orientation","horizontal","Location","southoutside");
box on
legend box off
fontsize(gcf,20,"points");
set(gcf,'position',[10,10,600,400])
saveas(gcf,'Absorption_LD2.png');


degree_D2=5/20;
rate_Absorption_D2=change_area./(4*degree_D2*area);


figure(10);
hold on
plot(f,rate_Absorption_D2,":","LineWidth",5);
hold off

rate_Absorption_D2_1=change_area/(4*degree_D2);
figure(11);
hold on
plot(f,rate_Absorption_D2_1,"LineWidth",5);
hold off
%thickenall_big
clear change_area


%% -------------- Plate_1 ----------------------
figure("Name","Absorption_p1");
[Plate1,Plate2]=set_default;
alpha=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha,"LineWidth",5);hold on;

Plate1.p1=25*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha);

Plate1.p1=20*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,":","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

Plate1.p1=35*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-.","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

Plate1.p1=40*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"--","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;


ylim([0,1.1]);
xticks(200:300:3000);  % Specify desired tick marks
set(gca, 'XTickLabel', arrayfun(@(v) sprintf('%g', v/1000), xticks, 'UniformOutput', false));
xlabel("Frequency in kHz");ylabel("Absorption coefficient");legend("3%","2.5%","2%","3.5%","4%","show","Orientation","horizontal","Location","southoutside");
box on
legend box off
fontsize(gcf,20,"points");
set(gcf,'position',[10,10,600,400])
saveas(gcf,'Absorption_p1.png');


degree_p1=0.005/0.03;
rate_Absorption_p1=change_area./(4*degree_p1*area);


figure(12);
hold on
plot(f,rate_Absorption_p1,"-","LineWidth",5);
hold off

rate_Absorption_p1_1=change_area/(4*degree_p1);
figure(13);
hold on
plot(f,rate_Absorption_p1_1,"LineWidth",5);
hold off
%thickenall_big
clear change_area



%%
figure("Name","Absorption_R1");
[Plate1,Plate2]=set_default;
alpha=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha,"LineWidth",5);hold on;

Plate1.R1=9*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha);

Plate1.R1=8*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,":","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

Plate1.R1=11*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-.","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

Plate1.R1=12*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"--","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;


ylim([0,1.1]);
xticks(200:300:3000);  % Specify desired tick marks
set(gca, 'XTickLabel', arrayfun(@(v) sprintf('%g', v/1000), xticks, 'UniformOutput', false));
xlabel("Frequency in kHz");ylabel("Absorption coefficient");legend("10mm","9mm","8mm","11mm","12mm","show","Orientation","horizontal","Location","southoutside");
box on
legend box off
fontsize(gcf,20,"points");
set(gcf,'position',[10,10,600,400]);
saveas(gcf,'Absorption_R1.png');



degree_R1=1/10;
rate_Absorption_R1=change_area./(4*degree_R1*area);


figure(12);
hold on
plot(f,rate_Absorption_R1,"LineWidth",5);
hold off

rate_Absorption_R1_1=change_area/(4*degree_R1);
figure(13);
hold on
plot(f,rate_Absorption_R1_1,"--","LineWidth",5);
hold off
%thickenall_big
clear change_area

%%
figure("Name","Absorption_D1");
[Plate1,Plate2]=set_default;
alpha=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha,"LineWidth",5);hold on;

Plate1.D1=25*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha);

Plate1.D1=20*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,":","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

Plate1.D1=35*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-.","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

Plate1.D1=40*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"--","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;


ylim([0,1.1]);
xticks(200:300:3000);  % Specify desired tick marks
set(gca, 'XTickLabel', arrayfun(@(v) sprintf('%g', v/1000), xticks, 'UniformOutput', false));
xlabel("Frequency in kHz");ylabel("Absorption coefficient");legend("30mm","25mm","20mm","35mm","40mm","show","Orientation","horizontal","Location","southoutside");
box on
legend box off
fontsize(gcf,20,"points");
set(gcf,'position',[10,10,600,400])
saveas(gcf,'Absorption_LD1.png');


degree_D1=5/30;
rate_Absorption_D1=change_area./(4*degree_D1*area);


figure(12);
hold on
plot(f,rate_Absorption_D1,"-.","LineWidth",5);
hold off

rate_Absorption_D1_1=change_area/(4*degree_D1);
figure(13);
hold on
plot(f,rate_Absorption_D1_1,"LineWidth",5);
hold off
%thickenall_big
clear change_area

%%
figure("Name","Absorption_d1");
[Plate1,Plate2]=set_default;
alpha=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha,"LineWidth",5);hold on;

Plate1.d1=0.4*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha);

Plate1.d1=0.3*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,":","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

Plate1.d1=0.6*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"-.","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;

Plate1.d1=0.7*10^(-3);
alpha_changed=Absorption_DLMMPps(f,envir_cons,Plate1,Plate2);
plot(f,alpha_changed,"--","LineWidth",5);hold on;
change_area=deltaf*abs(alpha_changed-alpha)+change_area;


ylim([0,1.1]);
xticks(200:300:3000);  % Specify desired tick marks
set(gca, 'XTickLabel', arrayfun(@(v) sprintf('%g', v/1000), xticks, 'UniformOutput', false));
xlabel("Frequency in kHz");ylabel("Absorption coefficient");legend("0.5mm","0.4mm","0.3mm","0.6mm","0.7mm","show","Orientation","horizontal","Location","southoutside");
box on
legend box off
fontsize(gcf,20,"points");
set(gcf,'position',[10,10,600,400])
saveas(gcf,'Absorption_d1.png');


degree_d1=0.1/0.5;
rate_Absorption_d1=change_area./(4*degree_d1*area);


figure(12);
hold on
plot(f,rate_Absorption_d1,":","LineWidth",5);
hold off

rate_Absorption_d1_1=change_area/(4*degree_d1);
figure(13);
hold on
plot(f,rate_Absorption_d1_1,"LineWidth",5);
hold off
%thickenall_big
clear change_area

%% legend

figure(10);
xlabel("frequency in Hz");ylabel("Sensitivity");
legend("t","R3","d3","d2","p2","p3","D2");
fontsize(gcf,20,"points");
saveas(gcf,'Sensitivity_Plate2.png');

figure(11);
xlabel("frequency in Hz");ylabel("Sensitivity");
legend("t","R3","d3","d2","p2","p3","D2");
fontsize(gcf,20,"points");


figure(12);
xlabel("frequency in Hz");ylabel("Sensitivity");
legend("p1","R1","D1","d1");
grid on
fontsize(gcf,20,"points");
saveas(gcf,'Sensitivity_Plate1.png');

figure(13);
xlabel("frequency in Hz");ylabel("Sensitivity");
legend("p1","R1","D1","d1");
grid on
fontsize(gcf,20,"points");

%% quantification
%sensitivity_Absorption_d1 = mean(rate_Absorption_d1(590:3277)); %500-2000HZ
sensitivity_Absorption_d1 = mean(rate_Absorption_d1(590:2415));% 500-1500HZ
sensitivity_Absorption_t = mean(rate_Absorption_t(590:2415));% 500-1500HZ
sensitivity_Absorption_R1 = mean(rate_Absorption_R1(590:2415));% 500-1500HZ
sensitivity_Absorption_p1 = mean(rate_Absorption_p1(590:2415));% 500-1500HZ
sensitivity_Absorption_D1 = mean(rate_Absorption_D1(590:2415));% 500-1500HZ

sensitivity_Absorption_d2 = mean(rate_Absorption_d2(590:2415));% 500-1500HZ
sensitivity_Absorption_d3 = mean(rate_Absorption_d3(590:2415));% 500-1500HZ
sensitivity_Absorption_R3 = mean(rate_Absorption_R3(590:2415));% 500-1500HZ
sensitivity_Absorption_p2 = mean(rate_Absorption_p2(590:2415));% 500-1500HZ
sensitivity_Absorption_p3 = mean(rate_Absorption_p3(590:2415));% 500-1500HZ
sensitivity_Absorption_D2 = mean(rate_Absorption_D2(590:2415));% 500-1500HZ


%% save figures
% FolderName = "/Users/elliot/Documents/Master thesis/Thereotical calculation/DLMMP/sweep";   % Your destination folder
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = get(FigHandle, 'Name');
%   savefig(FigHandle, fullfile(FolderName, FigName, '.fig'));
% end

%%
% sensitivity_Absorption.d1_2k = mean(rate_Absorption_d1(590:3277));% 500-1500HZ
% sensitivity_Absorption.t_2k  = mean(rate_Absorption_t(590:3277));% 500-1500HZ
% sensitivity_Absorption.R1_2k  = mean(rate_Absorption_R1(590:3277));% 500-1500HZ
% sensitivity_Absorption.p1_2k  = mean(rate_Absorption_p1(590:3277));% 500-1500HZ
% sensitivity_Absorption.D1_2k  = mean(rate_Absorption_D1(590:3277));% 500-1500HZ
% 
% sensitivity_Absorption.d2_2k  = mean(rate_Absorption_d2(590:3277));% 500-1500HZ
% sensitivity_Absorption.d3_2k  = mean(rate_Absorption_d3(590:3277));% 500-1500HZ
% sensitivity_Absorption.R3_2k  = mean(rate_Absorption_R3(590:3277));% 500-1500HZ
% sensitivity_Absorption.p2_2k  = mean(rate_Absorption_p2(590:3277));% 500-1500HZ
% sensitivity_Absorption.p3_2k  = mean(rate_Absorption_p3(590:3277));% 500-1500HZ
% sensitivity_Absorption.D2_2k  = mean(rate_Absorption_D2(590:3277));% 500-1500HZ

%%
function [Plate1,Plate2]=set_default()
    %hole diameter
    Plate1.d1=0.5*10^(-3);
    Plate2.d2=0.5*10^(-3);
    Plate2.d3=0.5*10^(-3);%holes in center
    %plate radius
    Plate1.R1=10*10^(-3);
    Plate2.R3=5*10^(-3);
    Plate2.R2=Plate1.R1-Plate2.R3;
    %thickness
    Plate1.t=1*10^(-3);
    %percentage of open area
    Plate1.p1=0.03;
    Plate2.p2=0.01;
    Plate2.p3=0.01;
    %depth in m 
    Plate1.D1=30*10^(-3);
    Plate2.D2=20*10^(-3);
end