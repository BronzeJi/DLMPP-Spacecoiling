function [alpha]=Absorption_DLMMP(f,envir_cons,Plate1,Plate2)
    w=2*pi*f;%omega
    c=envir_cons.c;
    u=envir_cons.mu;
    rho_c=envir_cons.rho_c;

    %hole diameter
    d1=Plate1.d1;
    d2=Plate2.d2;
    d3=Plate2.d3;
    %plate diameter
    R1=Plate1.R1;
    R3=Plate2.R3;
    %R2=R1-R3;
    %thickness
    t=Plate1.t;
    %percentage of open area
    p1=Plate1.p1;
    p2=Plate2.p2;
    p3=Plate2.p3;
    %depth in m 
    D1=Plate1.D1;
    D2=Plate2.D2;
    
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
    
    Z1=Z_mmp1+Z_c1;
    Z2=Z_mmp2+Z_c2;
    Z3=Z_mmp3+Z_c3;
    
    Z_1stlayer=Z1;
    Z_2ndlayer=(a2+a3)./((a3./Z3)+(a2./Z2));
    Z=Z_mmp1+((Z_c1.*Z_2ndlayer)./(Z_c1+Z_2ndlayer));
    
    
    %absortion coefficient
    Alpha=@(Za)(4*real(Za)./((1+real(Za)).^2+imag(Za).^2));
    absoropeZ1=Alpha(Z_1stlayer);% mass and air gap
%     absoropeZ2=Alpha(Z2);% mass and air gap
%     absoropeZ3=Alpha(Z3);% mass and air gap
    absoropeZ2ndlayer=Alpha(Z_2ndlayer);% mass and air gap
    absoropeZ=Alpha(Z);% mass and air gap
    
    
%     figure;
%     plot(f,absoropeZ1);hold on;
%     plot(f,absoropeZ2ndlayer); hold on;
%     plot(f,absoropeZ);hold off;
%     legend("First layer","Second layer","All","Location","best",'Location','eastoutside');
%     xlabel("Frequency in Hz");
%     ylabel("Absorption coefficient");
%     grid on
    alpha=absoropeZ;
end

