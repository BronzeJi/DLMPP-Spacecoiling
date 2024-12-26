function [alpha_t,alpha0]=space_coiling_unequalopening_CR(f,a,b1,h,n,d,CR)
    tic;
    %environment
    eta=1.8134*10^(-5);%viscosity
    P_0=1.013*10^5; % air pressure
    C_v=0.712;%kJ/(kg.K)
    gamma=1.4;%the ratio of specific heat Cp/Cv % CV and CP are the heat capacity at constant volume and pressure
    C_p=C_v*gamma;
    k_heat=0.026;%W/(m.k) thermal conductivity
    rho_0=1.21;%kg*m-3
    K_0=gamma*P_0;%adiabatic bulk modulus
    omega=2*pi*f;
    c_0=343;%m/s
    k_0=2*pi*f/c_0;
    Pr=C_p*eta/k_heat;%Pandtl number
    A0=(2*h+d)*a;
    A1=h*a;
    A00=A1;
    
    
    %different length 
    i=1:n;
%     b=b1*CR.^(1-i);
%     b=[b,b1*CR.^((i+n)-2*n)];

    b=b1.*CR.^(i-1);
    b=[b,b1*CR.^(n-i)];

    l=sqrt((h-b).^2+(b+d).^2);
    Length=sum(b);
    clear i;

    m_lim=0:10;
    n_lim=0:10;
    
    % Create matrices for i and j
    alpha=(2*m_lim+1)*pi/a;
    beta_unequal=(2*n_lim+1)*pi./b';
    
    G_rho=sqrt(1j*omega*rho_0/eta);
    G_K=sqrt(1j*omega*Pr*rho_0/eta);
    
    %K_f & rho_f
    
    K_part1=zeros(1,numel(f));
    rho_part1=zeros(1,numel(f));

    %-----------------------------------
    T = ones(2, 2, numel(f)); % Preallocate T
    for j=1:length(b)
        K_part2=(4*(gamma-1)*G_K.^2)./((a/2)^2*(b(j)/2)^2);
        beta=beta_unequal(j,:);
        for i=1:numel(f)
            K_part1(i)=sum(sum(1./((alpha.^2).*(beta.^2).*(alpha.^2+beta.^2+G_K(i)^2))));
            rho_part1(i)=sum(sum(1./((alpha.^2).*(beta.^2).*(alpha.^2+beta.^2+G_rho(i)^2))));
        end
        
        rho_f=rho_0*(a/2)^2*(b(j)/2)^2./(4*(G_rho.^2).*rho_part1);
        K_f=K_0./(gamma-K_part2.*K_part1);
        
        
        %Impedance
        Z=sqrt(K_f.*rho_f)/(a*b(j));% impedance divided by area
        k=omega./sqrt(K_f./rho_f);% bulk modulus/ wavenumber
    
        
        elapsed_time = toc;
        disp(['Elapsed time: ', num2str(elapsed_time), ' seconds']);
        tic;
        %Transfer matrix
        %T = ones(2, 2, numel(f)); % Preallocate T
        T_temp = zeros(2, 2, numel(f)); % Preallocate T_temp
        for i = 1:numel(f)
            T_temp(:, :, i) = [cos(k(i)*l(j)),1i*Z(i)*sin(k(i)*l(j));1i/Z(i)*sin(k(i)*l(j)),cos(k(i)*l(j))];
            if j==1
                T(:, :, i)=T_temp(:, :, i);
            elseif j<=length(b)
                T(:, :, i)=T(:, :, i)*T_temp(:, :, i);
            else
                disp("error section");
            end
        end

        if j==1
            disp("start");
                Z_1=Z;
        elseif j==length(b)
            disp("end");
                Z_2N=Z;
        end
    end
    %-----------------------------------


    %absorption coefficient
    %Alpha=@(Za)(4*real(Za)./((1+real(Za)).^2+imag(Za).^2));
    %alpha=Alpha(Z);
    
    % Stop the timer and display the elapsed time
    elapsed_time = toc;
    disp(['Elapsed time: ', num2str(elapsed_time), ' seconds']);
    T_12=reshape(T(1,2,:),[1,numel(f)]);
    T_11=reshape(T(1,1,:),[1,numel(f)]);
    T_21=reshape(T(2,1,:),[1,numel(f)]);
    T_22=reshape(T(2,2,:),[1,numel(f)]);
    
    %T=(2*Z_1)./(T_11.*Z_2N+T_12+T_21.*Z_1.*Z_2N+T_22.*Z_1);
    %T= T_22;
    %Matrix solution for absorption    
    %M = zeros(2, 2, numel(f)); % Preallocate T
    Pvector0=zeros(6,numel(f));
    Pvector=zeros(8,numel(f));
    for i = 1:numel(f)
        M=[0,1,1,(T_12(i)/Z_2N(i))-T_11(i),-(T_12(i)/Z_2N(i))-T_11(i),0,0,0;
            0,1,-1,Z_1(i)*(T_21(i)-T_22(i)/Z_2N(i)),Z_1(i)*(T_21(i)+T_22(i)/Z_2N(i)),0,0,0;
            -1,1,1,0,0,0,0,0;
            -1,0,0,0,0,1,1,0;
            1,-rho_0*c_0/(Z_1(i)*A0),+rho_0*c_0/(Z_1(i)*A0),0,0,-A1/A0,A1/A0,0;
            0,0,0,1,1,-exp(1i*k_0(i)*Length),-exp(-1i*k_0(i)*Length),0;
            0,0,0,0,0,exp(1i*k_0(i)*Length),exp(-1i*k_0(i)*Length),-1;
            0,0,0,-1/Z_2N(i),1/Z_2N(i),-(A1/rho_0*c_0)*exp(1i*k_0(i)*Length),(A1/rho_0*c_0)*exp(-1i*k_0(i)*Length),-(A0/rho_0*c_0)];
    
    
        M0=[0,1,1,-T_11(i)+(T_12(i)/Z_2N(i)),-T_11(i)-(T_12(i)/Z_2N(i)),0;...
           0,1,-1,Z_1(i)*(T_21(i)-T_22(i)/Z_2N(i)),Z_1(i)*(T_21(i)+T_22(i)/Z_2N(i)),0;...
           -1,1,1,0,0,0;...
           1,-rho_0*c_0/(Z_1(i)*A00),rho_0*c_0/(Z_1(i)*A00),0,0,0;...
           0,0,0,-1/Z_2N(i),1/Z_2N(i),-A00/(rho_0*c_0);...
           0,0,0,1,1,-1];


        Pvector0(:,i)=M0\[0;0;1;1;0;0];
        Pvector(:,i)=M\[0,0,1,1,1,0,0,0]';
    end
    alpha0=1-(abs(Pvector0(6,:)/1)).^2;
    alpha_t=1-(abs(Pvector(8,:)/1)).^2;
end