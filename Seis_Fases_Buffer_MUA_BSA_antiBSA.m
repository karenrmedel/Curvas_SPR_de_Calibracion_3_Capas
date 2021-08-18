clear all
format long

lambda=632.8e-9; %Longitud de onda (nm)

% DECLARACIÓN DE PARÁMETROS DE LAS MONOCAPAS
% Los índices de refracción (el número complejo n + k) y los grosores de
% las distintas capas (en nanómetros) son registrados en primer lugar.
% Para las capas de los extremos, sólo n y k.
% En este caso, se utilizan ocho materiales distintos para un sistema 
% máximo de seis capas.

% 1: Prisma de vidrio SF14
n_prisma=1.723;
k_prisma=0;

% 2: Oro
n_oro=0.1726; %Valor en 633 nm
k_oro=3.4218; %Valor en 633 nm
d_oro=47; %Valor arbitrario para una película delgada comercial (47 nm)

% 3: Monocapa 1 (MUA)
n_pr1=1.483; %Índice de refracción de MUA
k_pr1=0;
d_pr1=15.23*cos(pi/6); %Tamaño de MUA (15.23 A)
% El coseno es debido a la inclinación natural sobre una monocapa en oro.

% 4: Monocapa 2 (EDC)
n_pr5=1.459; %Índice de refracción de NHS
k_pr5=0;
d_pr5=2.104*cos(pi/6); %Tamaño de EDC (2.104 A)

% 5: Monocapa 2 (NHS)
n_pr4=1.6; %Índice de refracción de NHS
k_pr4=0;
d_pr4=5.76*cos(pi/6); %Tamaño de NHS(5.76 A)

% 6: Monocapa 2 (BSA)
n_pr2=1.575; %Índice de refracción de BSA
k_pr2=0;
d1_pr2_1=61.31; %(61.31 A Mas probable)
d_pr2_1=d1_pr2_1*cos(pi/6);

% 7: Monocapa 3 (antiBSA)
n_pr3=1.542; %Índice de refracción de antiBSA
k_pr3=0;
d_pr3=40; %Tamaño de antiBSA (40 A)

% 8: Buffer Acetato de Sodio (medio)
n_med=1.37;
k_med=0;

% VIDRIO/ORO/BUFFER

% Creación de matrices para n, k y d y cálculo de la constante
% dieléctrica de cada capa

en_C(1)=n_prisma;
ek_C(1)=0;

en_C(2)=n_oro;
ek_C(2)=k_oro;
d_1_C(2)=1e-9*(d_oro); % Grosor en nm

en_C(3)=n_med;
ek_C(3)=k_med;

for x=1:length(en_C);
    er_C=en_C(x)^2-ek_C(x)^2;
    ei_C=2*en_C(x)*ek_C(x);
    e_C(x)=complex(er_C,ei_C);
end

% *************************************************************************
% CÁLCULOS
theta_deg=65:.01:80; % Rango: 50° a 80°. Intervalos: 0.01°
theta_ext=theta_deg/180*pi;
theta0=(pi/4)+asin(1/en_C(1)*sin(theta_ext-(pi/4)));

for jtheta=1:length(theta0);
    theta=theta0(jtheta);
    q1=sqrt(e_C(1)-en_C(1)^2*sin(theta)^2)/e_C(1);
    qn=sqrt(e_C(end)-en_C(1)^2*sin(theta)^2)/e_C(end);
    
    for j=2:(length(e_C)-1)
        beta=d_1_C(j)*2*pi/lambda*sqrt(e_C(j)-en_C(1)^2*sin(theta)^2);
        q=sqrt(e_C(j)-en_C(1)^2*sin(theta)^2)/e_C(j);
        em(j,1,1)=cos(beta);
        em(j,1,2)=-i*sin(beta)/q;
        em(j,2,1)=-i*sin(beta)*q;
        em(j,2,2)=cos(beta);
    end
    
    emtot=[1 0; 0 1];
    
    for j=2:(length(e_C)-1)
        emtot1(:,:)=em(j,:,:);
        emtot=emtot*emtot1;
    end
    
    rp=((emtot(1,1)+emtot(1,2)*qn)*q1-(emtot(2,1)+emtot(2,2)*qn))/...
        ((emtot(1,1)+emtot(1,2)*qn)*q1+(emtot(2,1)+emtot(2,2)*qn));
    ref_1_C=rp*conj(rp);
    refle_1_C(jtheta)=ref_1_C*100;
end

refmin_1_C=find(refle_1_C==min(refle_1_C));
thetacr_1_C=theta_deg(refmin_1_C); % Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico es: ');
disp(thetacr_1_C);
disp('Para d =');
disp(d1_pr2_1);
% *************************************************************************


% VIDRIO/ORO/MUA/BUFFER

% Creación de matrices para n, k y d y cálculo de la constante
% dieléctrica de cada capa

en_A(1)=n_prisma;
ek_A(1)=0;

en_A(2)=n_oro;
ek_A(2)=k_oro;
d_1_A(2)=1e-9*(d_oro); % Grosor en nm

en_A(3)=n_pr1;
ek_A(3)=k_pr1;
d_1_A(3)=1e-10*(d_pr1); % Grosor en A

en_A(4)=n_med;
ek_A(4)=k_med;

for x=1:length(en_A);
    er_A=en_A(x)^2-ek_A(x)^2;
    ei_A=2*en_A(x)*ek_A(x);
    e_A(x)=complex(er_A,ei_A);
end

% *************************************************************************
% CÁLCULOS
theta0=(pi/4)+asin(1/en_A(1)*sin(theta_ext-(pi/4)));

for jtheta=1:length(theta0);
    theta=theta0(jtheta);
    q1=sqrt(e_A(1)-en_A(1)^2*sin(theta)^2)/e_A(1);
    qn=sqrt(e_A(end)-en_A(1)^2*sin(theta)^2)/e_A(end);
    
    for j=2:(length(e_A)-1)
        beta=d_1_A(j)*2*pi/lambda*sqrt(e_A(j)-en_A(1)^2*sin(theta)^2);
        q=sqrt(e_A(j)-en_A(1)^2*sin(theta)^2)/e_A(j);
        em(j,1,1)=cos(beta);
        em(j,1,2)=-i*sin(beta)/q;
        em(j,2,1)=-i*sin(beta)*q;
        em(j,2,2)=cos(beta);
    end
    
    emtot=[1 0; 0 1];
    
    for j=2:(length(e_A)-1)
        emtot1(:,:)=em(j,:,:);
        emtot=emtot*emtot1;
    end
    
    rp=((emtot(1,1)+emtot(1,2)*qn)*q1-(emtot(2,1)+emtot(2,2)*qn))/...
        ((emtot(1,1)+emtot(1,2)*qn)*q1+(emtot(2,1)+emtot(2,2)*qn));
    ref_1_A=rp*conj(rp);
    refle_1_A(jtheta)=ref_1_A*100;
end

refmin_1_A=find(refle_1_A==min(refle_1_A));
thetacr_1_A=theta_deg(refmin_1_A); % Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico es: ');
disp(thetacr_1_A);
disp('Para d =');
disp(d1_pr2_1);
% *************************************************************************


% VIDRIO/ORO/MUA/EDC/BUFFER

% Creación de matrices para n, k y d y cálculo de la constante
% dieléctrica de cada capa

en_D(1)=n_prisma;
ek_D(1)=0;

en_D(2)=n_oro;
ek_D(2)=k_oro;
d_1_D(2)=1e-9*(d_oro); % Grosor en nm Oro

en_D(3)=n_pr1;
ek_D(3)=k_pr1;
d_1_D(3)=1e-10*(d_pr1); % Grosor en A MUA

en_D(4)=n_pr5;
ek_D(4)=k_pr5;
d_1_D(4)=1e-10*(d_pr5); % Grosor en A EDC

en_D(5)=n_med;
ek_D(5)=k_med;

for x=1:length(en_D);
    er_D=en_D(x)^2-ek_D(x)^2;
    ei_D=2*en_D(x)*ek_D(x);
    e_D(x)=complex(er_D,ei_D);
end

% *************************************************************************
% CÁLCULOS
theta0=(pi/4)+asin(1/en_D(1)*sin(theta_ext-(pi/4)));

for jtheta=1:length(theta0);
    theta=theta0(jtheta);
    q1=sqrt(e_D(1)-en_D(1)^2*sin(theta)^2)/e_D(1);
    qn=sqrt(e_D(end)-en_D(1)^2*sin(theta)^2)/e_D(end);
    
    for j=2:(length(e_D)-1)
        beta=d_1_D(j)*2*pi/lambda*sqrt(e_D(j)-en_D(1)^2*sin(theta)^2);
        q=sqrt(e_D(j)-en_D(1)^2*sin(theta)^2)/e_D(j);
        em(j,1,1)=cos(beta);
        em(j,1,2)=-i*sin(beta)/q;
        em(j,2,1)=-i*sin(beta)*q;
        em(j,2,2)=cos(beta);
    end
    
    emtot=[1 0; 0 1];
    
    for j=2:(length(e_D)-1)
        emtot1(:,:)=em(j,:,:);
        emtot=emtot*emtot1;
    end
    
    rp=((emtot(1,1)+emtot(1,2)*qn)*q1-(emtot(2,1)+emtot(2,2)*qn))/...
        ((emtot(1,1)+emtot(1,2)*qn)*q1+(emtot(2,1)+emtot(2,2)*qn));
    ref_1_D=rp*conj(rp);
    refle_1_D(jtheta)=ref_1_D*100;
end

refmin_1_D=find(refle_1_D==min(refle_1_D));
thetacr_1_D=theta_deg(refmin_1_D); % Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico es: ');
disp(thetacr_1_D);
disp('Para d =');
disp(d_pr5);
% *************************************************************************

% VIDRIO/ORO/MUA/NHS/BUFFER
% Creación de matrices para n, k y d y cálculo de la constante
% dieléctrica de cada capa

en_E(1)=n_prisma;
ek_E(1)=0;

en_E(2)=n_oro;
ek_E(2)=k_oro;
d_1_E(2)=1e-9*(d_oro); % Grosor en nm Oro

en_E(3)=n_pr1;
ek_E(3)=k_pr1;
d_1_E(3)=1e-10*(d_pr1); % Grosor en A MUA

en_E(4)=n_pr4;
ek_E(4)=k_pr4;
d_1_E(4)=1e-10*(d_pr4); % Grosor en A NHS
en_E(5)=n_med;
ek_E(5)=k_med;

for x=1:length(en_E);
    er_E=en_E(x)^2-ek_E(x)^2;
    ei_E=2*en_E(x)*ek_E(x);
    e_E(x)=complex(er_E,ei_E);
end

% *************************************************************************
% CÁLCULOS
theta0=(pi/4)+asin(1/en_E(1)*sin(theta_ext-(pi/4)));

for jtheta=1:length(theta0);
    theta=theta0(jtheta);
    q1=sqrt(e_E(1)-en_E(1)^2*sin(theta)^2)/e_E(1);
    qn=sqrt(e_E(end)-en_E(1)^2*sin(theta)^2)/e_E(end);
    
    for j=2:(length(e_E)-1)
        beta=d_1_E(j)*2*pi/lambda*sqrt(e_E(j)-en_E(1)^2*sin(theta)^2);
        q=sqrt(e_E(j)-en_E(1)^2*sin(theta)^2)/e_E(j);
        em(j,1,1)=cos(beta);
        em(j,1,2)=-i*sin(beta)/q;
        em(j,2,1)=-i*sin(beta)*q;
        em(j,2,2)=cos(beta);
    end
    
    emtot=[1 0; 0 1];
    
    for j=2:(length(e_E)-1)
        emtot1(:,:)=em(j,:,:);
        emtot=emtot*emtot1;
    end
    
    rp=((emtot(1,1)+emtot(1,2)*qn)*q1-(emtot(2,1)+emtot(2,2)*qn))/...
        ((emtot(1,1)+emtot(1,2)*qn)*q1+(emtot(2,1)+emtot(2,2)*qn));
    ref_1_E=rp*conj(rp);
    refle_1_E(jtheta)=ref_1_E*100;
end

refmin_1_E=find(refle_1_E==min(refle_1_E));
thetacr_1_E=theta_deg(refmin_1_E); % Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico es: ');
disp(thetacr_1_E);
disp('Para d =');
disp(d_pr4);
% *************************************************************************

% VIDRIO/ORO/MUA/BSA/BUFFER
% Creación de matrices para n, k y d y cálculo de la constante
% dieléctrica de cada capa

en_B(1)=n_prisma;
ek_B(1)=0;

en_B(2)=n_oro;
ek_B(2)=k_oro;
d_1_B(2)=1e-9*(d_oro); % Grosor en nm

en_B(3)=n_pr1;
ek_B(3)=k_pr1;
d_1_B(3)=1e-10*(d_pr1); % Grosor en A

en_B(4)=n_pr2;
ek_B(4)=k_pr2;
d_1_B(4)=1e-10*(d_pr2_1); % Grosor en A

en_B(5)=n_med;
ek_B(5)=k_med;

for x=1:length(en_B);
    er_B=en_B(x)^2-ek_B(x)^2;
    ei_B=2*en_B(x)*ek_B(x);
    e_B(x)=complex(er_B,ei_B);
end

% *************************************************************************
% CÁLCULOS
theta0=(pi/4)+asin(1/en_B(1)*sin(theta_ext-(pi/4)));

for jtheta=1:length(theta0);
    theta=theta0(jtheta);
    q1=sqrt(e_B(1)-en_B(1)^2*sin(theta)^2)/e_B(1);
    qn=sqrt(e_B(end)-en_B(1)^2*sin(theta)^2)/e_B(end);
    
    for j=2:(length(e_B)-1)
        beta=d_1_B(j)*2*pi/lambda*sqrt(e_B(j)-en_B(1)^2*sin(theta)^2);
        q=sqrt(e_B(j)-en_B(1)^2*sin(theta)^2)/e_B(j);
        em(j,1,1)=cos(beta);
        em(j,1,2)=-i*sin(beta)/q;
        em(j,2,1)=-i*sin(beta)*q;
        em(j,2,2)=cos(beta);
    end
    
    emtot=[1 0; 0 1];
    
    for j=2:(length(e_B)-1)
        emtot1(:,:)=em(j,:,:);
        emtot=emtot*emtot1;
    end
    
    rp=((emtot(1,1)+emtot(1,2)*qn)*q1-(emtot(2,1)+emtot(2,2)*qn))/...
        ((emtot(1,1)+emtot(1,2)*qn)*q1+(emtot(2,1)+emtot(2,2)*qn));
    ref_1_B=rp*conj(rp);
    refle_1_B(jtheta)=ref_1_B*100;
end

refmin_1_B=find(refle_1_B==min(refle_1_B));
thetacr_1_B=theta_deg(refmin_1_B); % Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico es: ');
disp(thetacr_1_B);
disp('Para d =');
disp(d1_pr2_1);
% *************************************************************************


% VIDRIO/ORO/MUA/BSA/antiBSA/BUFFER
% Creación de matrices para n, k y d y cálculo de la constante
% dieléctrica de cada capa

en(1)=n_prisma;
ek(1)=0;

en(2)=n_oro;
ek(2)=k_oro;
d_1(2)=1e-9*(d_oro); % Grosor en nm

en(3)=n_pr1;
ek(3)=k_pr1;
d_1(3)=1e-10*(d_pr1); % Grosor en A

en(4)=n_pr2;
ek(4)=k_pr2;
d_1(4)=1e-10*(d_pr2_1); % Grosor en A

en(5)=n_pr3;
ek(5)=k_pr3;
d_1(5)=1e-10*(d_pr3); % Grosor en A

en(6)=n_med;
ek(6)=k_med;

for x=1:length(en);
    er=en(x)^2-ek(x)^2;
    ei=2*en(x)*ek(x);
    e(x)=complex(er,ei);
end

% *************************************************************************
% CÁLCULOS
theta0=(pi/4)+asin(1/en(1)*sin(theta_ext-(pi/4)));

for jtheta=1:length(theta0);
    theta=theta0(jtheta);
    q1=sqrt(e(1)-en(1)^2*sin(theta)^2)/e(1);
    qn=sqrt(e(end)-en(1)^2*sin(theta)^2)/e(end);
    
    for j=2:(length(e)-1)
        beta=d_1(j)*2*pi/lambda*sqrt(e(j)-en(1)^2*sin(theta)^2);
        q=sqrt(e(j)-en(1)^2*sin(theta)^2)/e(j);
        em(j,1,1)=cos(beta);
        em(j,1,2)=-i*sin(beta)/q;
        em(j,2,1)=-i*sin(beta)*q;
        em(j,2,2)=cos(beta);
    end
    
    emtot=[1 0; 0 1];
    
    for j=2:(length(e)-1)
        emtot1(:,:)=em(j,:,:);
        emtot=emtot*emtot1;
    end
    
    rp=((emtot(1,1)+emtot(1,2)*qn)*q1-(emtot(2,1)+emtot(2,2)*qn))/...
        ((emtot(1,1)+emtot(1,2)*qn)*q1+(emtot(2,1)+emtot(2,2)*qn));
    ref_1=rp*conj(rp);
    refle_1(jtheta)=ref_1*100;
end


refmin_1=find(refle_1==min(refle_1));
thetacr_1=theta_deg(refmin_1); % Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico es: ');
disp(thetacr_1);
disp('Para d =');
disp(d1_pr2_1);
% *************************************************************************

% DIBUJO DE LA GRÁFICA
figure(4)
clf
plot(theta_deg,refle_1_C,'-k', 'DisplayName',...
    ['Buffer,       \theta_{min} = ', num2str(thetacr_1_C), '°'])  
hold on;
plot(theta_deg,refle_1_A,':k', 'DisplayName',...
    ['MUA,         \theta_{min} = ', num2str(thetacr_1_A), '°']) 
plot(theta_deg,refle_1_D,'-c', 'DisplayName',...
    ['EDC,         \theta_{min} = ', num2str(thetacr_1_D), '°']) 
plot(theta_deg,refle_1_E,'--c', 'DisplayName',...
    ['NHS,         \theta_{min} = ', num2str(thetacr_1_E), '°'])
plot(theta_deg,refle_1_B,'-.k', 'DisplayName',...
    ['BSA,         \theta_{min} = ', num2str(thetacr_1_B), '°']) 
plot(theta_deg,refle_1,'--k', 'DisplayName',...
    ['antiBSA,   \theta_{min} = ', num2str(thetacr_1), '°']) 
hold off;

title('SIMULACION CURVA SPR MUA+BSA(61.31 A)+antiBSA(40 A)+Acetato de Potasio')
ylabel('Reflectancia (%)')
xlabel('Ángulo (°)')