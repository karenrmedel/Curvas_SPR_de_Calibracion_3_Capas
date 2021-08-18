clear all
format long

lambda=632.8e-9; %Longitud de onda (nm)

% DECLARACIÓN DE PARÁMETROS DE LAS MONOCAPAS
% Los índices de refracción (el número complejo n + k) y los grosores de
% las distintas capas (en nanómetros) son registrados en primer lugar.
% Para las capas de los extremos, sólo n y k.
% En este caso, se utilizan cinco materiales distintos.

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

% 4: Monocapa 2 (IgG)
n_pr2=1.542;
k_pr2=0;
d1_pr2_1=76.51;
d1_pr2_2=87.19;
d1_pr2_3=124.96;
d1_pr2_4=130.11;
d1_pr2_5=131.11;
d1_pr2_6=136.75;
d1_pr2_7=142.63; %(118.47 A Promedio)
d_pr2_1=d1_pr2_1*cos(pi/6);
d_pr2_2=d1_pr2_2*cos(pi/6);
d_pr2_3=d1_pr2_3*cos(pi/6);
d_pr2_4=d1_pr2_4*cos(pi/6);
d_pr2_5=d1_pr2_5*cos(pi/6);
d_pr2_6=d1_pr2_6*cos(pi/6);
d_pr2_7=d1_pr2_7*cos(pi/6);

% 5: Buffer HEPES (medio)
n_med=1.3341;
k_med=0;

% Creación de matrices para n, k y d y cálculo de la constante
% dieléctrica de cada capa

en(1)=n_prisma;
ek(1)=0;

en(2)=n_oro;
ek(2)=k_oro;
d_1(2)=1e-9*(d_oro); % Grosor en nm
d_2(2)=1e-9*(d_oro);
d_3(2)=1e-9*(d_oro);
d_4(2)=1e-9*(d_oro);
d_5(2)=1e-9*(d_oro);
d_6(2)=1e-9*(d_oro);
d_7(2)=1e-9*(d_oro);

en(3)=n_pr1;
ek(3)=k_pr1;
d_1(3)=1e-10*(d_pr1);
d_2(3)=1e-10*(d_pr1); % Grosor en nm
d_3(3)=1e-10*(d_pr1);
d_4(3)=1e-10*(d_pr1);
d_5(3)=1e-10*(d_pr1);
d_6(3)=1e-10*(d_pr1);
d_7(3)=1e-10*(d_pr1);

en(4)=n_pr2;
ek(4)=k_pr2;
d_1(4)=1e-10*(d_pr2_1);
d_2(4)=1e-10*(d_pr2_2); % Grosor en nm
d_3(4)=1e-10*(d_pr2_3);
d_4(4)=1e-10*(d_pr2_4);
d_5(4)=1e-10*(d_pr2_5);
d_6(4)=1e-10*(d_pr2_6);
d_7(4)=1e-10*(d_pr2_7);

en(5)=n_med;
ek(5)=k_med;

for x=1:length(en);
    er=en(x)^2-ek(x)^2;
    ei=2*en(x)*ek(x);
    e(x)=complex(er,ei);
end

% *************************************************************************
% CÁLCULOS
theta_deg=66:.01:75; % Rango: 50° a 80°. Intervalos: 0.01°
theta_ext=theta_deg/180*pi;
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

for jtheta=1:length(theta0);
    theta=theta0(jtheta);
    q1=sqrt(e(1)-en(1)^2*sin(theta)^2)/e(1);
    qn=sqrt(e(end)-en(1)^2*sin(theta)^2)/e(end);
    
    for j=2:(length(e)-1)
        beta=d_2(j)*2*pi/lambda*sqrt(e(j)-en(1)^2*sin(theta)^2);
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
    ref_2=rp*conj(rp);
    refle_2(jtheta)=ref_2*100;
end

for jtheta=1:length(theta0);
    theta=theta0(jtheta);
    q1=sqrt(e(1)-en(1)^2*sin(theta)^2)/e(1);
    qn=sqrt(e(end)-en(1)^2*sin(theta)^2)/e(end);
    
    for j=2:(length(e)-1)
        beta=d_3(j)*2*pi/lambda*sqrt(e(j)-en(1)^2*sin(theta)^2);
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
    ref_3=rp*conj(rp);
    refle_3(jtheta)=ref_3*100;
end

for jtheta=1:length(theta0);
    theta=theta0(jtheta);
    q1=sqrt(e(1)-en(1)^2*sin(theta)^2)/e(1);
    qn=sqrt(e(end)-en(1)^2*sin(theta)^2)/e(end);
    
    for j=2:(length(e)-1)
        beta=d_4(j)*2*pi/lambda*sqrt(e(j)-en(1)^2*sin(theta)^2);
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
    ref_4=rp*conj(rp);
    refle_4(jtheta)=ref_4*100;
end

for jtheta=1:length(theta0);
    theta=theta0(jtheta);
    q1=sqrt(e(1)-en(1)^2*sin(theta)^2)/e(1);
    qn=sqrt(e(end)-en(1)^2*sin(theta)^2)/e(end);
    
    for j=2:(length(e)-1)
        beta=d_5(j)*2*pi/lambda*sqrt(e(j)-en(1)^2*sin(theta)^2);
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
    ref_5=rp*conj(rp);
    refle_5(jtheta)=ref_5*100;
end

for jtheta=1:length(theta0);
    theta=theta0(jtheta);
    q1=sqrt(e(1)-en(1)^2*sin(theta)^2)/e(1);
    qn=sqrt(e(end)-en(1)^2*sin(theta)^2)/e(end);
    
    for j=2:(length(e)-1)
        beta=d_6(j)*2*pi/lambda*sqrt(e(j)-en(1)^2*sin(theta)^2);
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
    ref_6=rp*conj(rp);
    refle_6(jtheta)=ref_6*100;
end

for jtheta=1:length(theta0);
    theta=theta0(jtheta);
    q1=sqrt(e(1)-en(1)^2*sin(theta)^2)/e(1);
    qn=sqrt(e(end)-en(1)^2*sin(theta)^2)/e(end);
    
    for j=2:(length(e)-1)
        beta=d_7(j)*2*pi/lambda*sqrt(e(j)-en(1)^2*sin(theta)^2);
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
    ref_7=rp*conj(rp);
    refle_7(jtheta)=ref_7*100;
end

refmin_1=find(refle_1==min(refle_1));
thetacr_1=theta_deg(refmin_1); % Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico es: ');
disp(thetacr_1);
disp('Para d =');
disp(d1_pr2_1);
% *************************************************************************

refmin_2=find(refle_2==min(refle_2));
thetacr_2=theta_deg(refmin_2); % Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico es: ');
disp(thetacr_2);
disp('Para d =');
disp(d1_pr2_2);
% *************************************************************************

refmin_3=find(refle_3==min(refle_3));
thetacr_3=theta_deg(refmin_3); % Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico es: ');
disp(thetacr_3);
disp('Para d =');
disp(d1_pr2_3);
% *************************************************************************

refmin_4=find(refle_4==min(refle_4));
thetacr_4=theta_deg(refmin_4); % Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico es: ');
disp(thetacr_4);
disp('Para d =');
disp(d1_pr2_4);
% *************************************************************************

refmin_5=find(refle_5==min(refle_5));
thetacr_5=theta_deg(refmin_5); % Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico es: ');
disp(thetacr_5);
disp('Para d =');
disp(d1_pr2_5);
% *************************************************************************

refmin_6=find(refle_6==min(refle_6));
thetacr_6=theta_deg(refmin_6); % Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico es: ');
disp(thetacr_6);
disp('Para d =');
disp(d1_pr2_6);
% *************************************************************************

refmin_7=find(refle_7==min(refle_7));
thetacr_7=theta_deg(refmin_7); % Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico es: ');
disp(thetacr_7);
disp('Para d =');
disp(d1_pr2_7);
% *************************************************************************

% DIBUJO DE LA GRÁFICA
figure(1)
clf
plot(theta_deg,refle_1, 'DisplayName',...
    ['d_{IgG} = ', num2str(d1_pr2_1), ' A, \theta_{min} = ', num2str(thetacr_1), '°'])  
hold on;
plot(theta_deg,refle_2, 'DisplayName',...
    ['d_{IgG} = ', num2str(d1_pr2_2), ' A, \theta_{min} = ', num2str(thetacr_2), '°']) 
plot(theta_deg,refle_3, 'DisplayName',...
    ['d_{IgG} = ', num2str(d1_pr2_3), ' A, \theta_{min} = ', num2str(thetacr_3), '°'])
plot(theta_deg,refle_4, 'DisplayName',...
    ['d_{IgG} = ', num2str(d1_pr2_4), ' A, \theta_{min} = ', num2str(thetacr_4), '°'])
plot(theta_deg,refle_5, 'DisplayName',...
    ['d_{IgG} = ', num2str(d1_pr2_5), ' A, \theta_{min} = ', num2str(thetacr_5), '°'])
plot(theta_deg,refle_6,'c', 'DisplayName',...
    ['d_{IgG} = ', num2str(d1_pr2_6), ' A, \theta_{min} = ', num2str(thetacr_6), '°'])
plot(theta_deg,refle_7,'m', 'DisplayName',...
    ['d_{IgG} = ', num2str(d1_pr2_7), ' A, \theta_{min} = ', num2str(thetacr_7), '°'])
hold off;

title('SIMULACION CURVA SPR MUA+IgG+HEPES')
ylabel('Reflectancia (%)')
xlabel('Ángulo (°)')