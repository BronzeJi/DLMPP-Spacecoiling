%% equal opening
a=50*10^(-3);
b=15*10^(-3);
h=52*10^(-3);
n=12;


d=1*10^(-3);%thickness
l=sqrt((h-b)^2+(b+d)^2);

f=linspace(10,3000,30000);
[alpha_t,alpha0]=space_coiling_equalopening(f,a,b,h,l,n,d);
harmo=zeros(40,1);
for i=1:40
    harmo(i)=i*271;
end

figure(2);
hold on
plot(f,alpha_t,"LineWidth",5);hold on;
%plot(f,alpha0,"--","LineWidth",4);hold off;
% xticks(200:200:2000);  % Specify desired tick marks
% set(gca, 'XTickLabel', arrayfun(@(v) sprintf('%g', v/1000), xticks, 'UniformOutput', false));
ylim([-0.3,1]);xlim([100,1900]);ylabel("Absorption coefficient");xlabel("Frequency in Hz");
%legend("b=35mm h=107 mm n=6","b=35mm h=83mm n=8","b=21mm h=64mm n=10","b=15mm h=52mm n=12mm","Location","southoutside","NumColumns", 2);
% grid on
% hold on; xline(harmo,"linewidth",4)
fontsize(gcf,20,"points");
% hold off
%thickenall_big

%% Unequal opening
clear
%close all
% %———————————————————————————— A perfect setting for Fan 1 60%
% a=100*10^(-3);
% b=30*10^(-3);
% h=90*10^(-3);
% n=9;
% %———————————————————————————— 1
a=50*10^(-3);
b=15*10^(-3);
h=40*10^(-3);
n=6;
%------------------------ 4
% a=50*10^(-3);
% b=10*10^(-3);
% h=55*10^(-3);
% n=10;
%------------------------ 2
% a=50*10^(-3);
% b=20*10^(-3)*1.5;
% h=40*10^(-3)*2;
% n=8;
% %-------------------3 with cr=1.5
% a=50*10^(-3);
% b=30*10^(-3);
% h=70*10^(-3);
% n=4;
d=1*10^(-3);%thickness
f=linspace(20,5000,30000);
CR=1.2;
[alpha_t,alpha0]=space_coiling_unequalopening_CR(f,a,b,h,n,d,CR);

% fan tonal noise
harmo=zeros(40,1);
for i=1:40
    harmo(i)=i*271;
end

figure(25);
hold on
plot(f,alpha_t,"LineWidth",5);hold on;
xline(harmo,"--","linewidth",3);
%plot(f,alpha0,"LineWidth");
hold off;
ylim([-0.1,1]);xlim([10,2000]);ylabel("Absorption coefficient");xlabel("Frequency in Hz");%legend("Equal openning cell with ventilation channels","Equal openning cell");
grid on
%legend("Unequal open cell with ventilation channels","Unequal open cell");

fontsize(gcf,20,"points");

%xlim([200,1500]);
%hold off
%% optimization
% TIME ISSUE

% clear
% close all
% % a=60*10^(-3);
% % b=15*10^(-3);
% % h=40*10^(-3);
% % n=6;
% % d=1*10^(-3);%thickness
% 
% nvars=3;
% lb=[20*10^(-3),30*10^(-3),5];% lower limit of parameter. frequency; top plate radius; neck length
% ub=[40*10^(-3),50*10^(-3),8];
% 
% tic;
% %options = optimoptions('particleswarm', 'MaxIterations', 1000);
% options = optimoptions('particleswarm', 'MaxTime', 60,'Display','iter','DisplayInterval',3); 
% [alpha_opt, fval] = particleswarm(@Optimization_Fano, nvars, lb, ub, options);
% disp(['Elapsed time total: ', num2str(elapsed_time), ' seconds']);
% % Display the optimal value of alpha
% disp(['Optimal parameter: ', num2str(alpha_opt)]);
% disp(['Objective function value: ', num2str(fval)]);

