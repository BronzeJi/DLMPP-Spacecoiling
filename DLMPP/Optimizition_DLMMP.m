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

%Frequency
f=linspace(100,3000,5000);

%Hole diameter
Plate1.d1=0.5*10^(-3);
Plate2.d2=0.5*10^(-3);
Plate2.d3=0.4*10^(-3);%holes in center
%Plate radius
Plate1.R1=12.5*10^(-3);%r2 in paper
Plate2.R3=9.5*10^(-3);
Plate2.R2=Plate1.R1-Plate2.R3;
%Thickness
Plate1.t=1*10^(-3);
%Prcentage of open area
Plate1.p1=0.03;
Plate2.p2=0.01;
Plate2.p3=0.01;
%Depth in m 
Plate1.D1=30*10^(-3);
Plate2.D2=20*10^(-3);

%Parameters to be optimized 
nvars2=4;
lb2=[3*10^(-3),0.2*10^(-3),0.009, 10*10^(-3)];% lower limit of parameter. frequency; top plate radius; neck length
ub2=[20*10^(-3),1*10^(-3),0.05, 30*10^(-3)];

options = optimoptions('particleswarm', 'MaxIterations', 1000);
[alpha_opt, fval] = particleswarm(@optimize, nvars2, lb2, ub2, options);%<<-------------
% Display the optimal value of alpha
disp(['Optimal parameter: ', num2str(alpha_opt)]);
disp(['Objective function value: ', num2str(fval)]);

Plate2.R3=alpha_opt(1);%<<-------------
Plate1.t=alpha_opt(2);%<<-------------
Plate1.p1=alpha_opt(3);%<<-------------
Plate2.D2=alpha_opt(4);%<<-------------


[absoropeZ,Z]=Absorption_DLMMP(f,envir_cons,Plate1,Plate2);


%%
% change the "f" to change the target frequency band
function [alpha]=optimize(x)% R3, t, D2
%optimize with a certain frequency band
    f=linspace(100,3000,10000);%<<-------------

    R3=x(1);
    t=x(2);%mm
    p1=x(3);%mm%<<-------------
    D2=x(4);

    Alpha=@(Za)(4*real(Za)./((1+real(Za)).^2+imag(Za).^2));
    w=2*pi*f;%omega  
    %environment constant
    c=343;%speed of sound in m/s
    u=1.84*10^(-5);%mue pa*s
    rho_c=1.213;% air density
    
    %hole diameter
    d1=0.5*10^(-3);%mm
    d2=0.5*10^(-3);%mm
    d3=0.5*10^(-3);
    %plate diameter
    R1=12.5*10^(-3);%r2 in paper
    %R3=9.5*10^(-3);
    %thickness
    %t=1*10^(-3);
    %percentage of open area
    %p1=0.03;
    p2=0.01;
    p3=0.01;
    %depth in m 
    D1=30*10^(-3);
    %D2=20*10^(-3);
    
    %-------------------------
    %area
    a1=R1^2*pi;
    a3=R3^2*pi;
    a2=a1-a3;
    
    x1=d1/2*sqrt(w/u);%x1=d*sqrt(f/10);viscosity parameter
    x2=d2/2*sqrt(w/u);%x1=d*sqrt(f/10);viscosity parameter
    x3=d3/2*sqrt(w/u);%x1=d*sqrt(f/10);viscosity parameter
    
    mass=@(p,x,d)(t/(p*c))*(1+1./sqrt(9+((x.^2)/2))+0.85*(d/t));
    resistor=@(p,x,d)32*u*t/(p*c*rho_c*d^2)*(sqrt(1+(x.^2)/32)+(x*d*sqrt(2)/(8*t)));
    cavity=@(D)-1i*cot(w*D/c);
    % equivalent mass and resistor accroding to Maa
    m1=mass(p1,x1,d1);
    r1=resistor(p1,x1,d1);
    Z_mmp1=r1+1j.*w.*m1;
    Z_c1=cavity(D1);
    
    % second layer. Be careful with the equivalent cavity depths! D=V/a
    V3=frustumConeVolume(R1,R3,D2);
    V2=D2*a1-V3;
    
    D_2=V2/a2;
    m2=mass(p2,x2,d2);
    r2=resistor(p2,x2,d2);
    Z_mmp2=r2+1j.*w.*m2;
    Z_c2=cavity(D_2);
    
    D_3=V3/a3;
    m3=mass(p3,x3,d3);
    r3=resistor(p3,x3,d3);
    Z_mmp3=r3+1j.*w.*m3;
    Z_c3=cavity(D_3);
    
    %Z1=Z_mmp1+Z_c1;
    Z2=Z_mmp2+Z_c2;
    Z3=Z_mmp3+Z_c3;
    
    Z_2ndlayer=(a2+a3)./((a3./Z3)+(a2./Z2));
    Z=Z_mmp1+((Z_c1.*Z_2ndlayer)./(Z_c1+Z_2ndlayer));
    
    %absortion coefficient
    absoropeZ=Alpha(Z);% mass and air gap
    alpha=mean(1-absoropeZ);

end