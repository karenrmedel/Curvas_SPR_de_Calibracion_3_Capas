clear all;
figure(1)
clf
clc 

%---PARÁMETROS DE LAS DISTINTAS CAPAS
%Vidrio (SF10), 633 nm
enpr=1.76;          %n del vidrio (sólo hay parte real)

%Oro, 633 nm
e1n=0.1726;         %Parte real de n (oro)
e1k=3.4218;         %Parte imaginaria de n (oro)
e1r=e1n^2-e1k^2;    %Parte real de epsilon (oro)
e1i=2*e1n*e1k;      %Parte imaginaria de epsilon (oro)
e1=complex(e1r,e1i); %Epsilon compleja (oro)
d1=47e-9;           %Grosor de la capa de oro (m)

%Agua
n2=1.333;               %n del analito (sólo hay parte real)
e2=n2^2;            %epsilon del analito (sólo hay parte real)

%---CONSTANTES Y PARÁMETROS DADOS 
c=2.99792458e8;
lambda=632.8e-9;
omega=2*pi/lambda*c;

%---ESCANEO DE LOS ÁNGULOS Y SOLUCIÓN DE LAS ECUACIONES DE FRESNEL
ang0=50;            %Límite inferior
ang1=62;            %Límite superior
% vals=10000;
% interval=ang1-ang0;
% angmat=ang0:(interval/vals):ang1; %Vector de ángulos
theta_deg=ang0:.01:ang1; % Rango: 48° a 60°. Intervalos: 0.01°  MATRIZ DE GRADOS DESDE 48 A 60 ( DE .01 EN .01)
theta_ext=theta_deg/180*pi; %Transformamos grados a radianes

for x=1:(length(theta_ext))
    theta=theta_ext(x);
	rpr1=(cos(theta)/enpr - sqrt(e1-(enpr^2)*(sin(theta)^2))/e1)/...
		(cos(theta)/enpr + sqrt(e1-(enpr^2)*(sin(theta)^2))/e1);
	r12=(sqrt(e1-(enpr^2)*(sin(theta)^2))/e1 -...
		sqrt(e2-(enpr^2)*(sin(theta)^2))/e2) /...
		(sqrt(e1-(enpr^2)*(sin(theta)^2))/e1 +...
		sqrt(e2-(enpr^2)*(sin(theta)^2))/e2);
	kz1d1=omega/c*d1*sqrt(e1-(enpr^2)*(sin(theta)^2));
	rpr12=(rpr1+ r12*exp(2*i*kz1d1))/(1+rpr1*r12*exp(2*i*kz1d1));
	ref(x)=rpr12*(conj(rpr12))*100;
end
 
refmin=find(ref==min(ref));
thetacr=theta_deg(refmin); %Posición del ángulo crítico (mín. reflectancia)

%
% SEGUNDA CURVA
%
%---PARÁMETROS DE LAS DISTINTAS CAPAS
%Vidrio (SF10), 633 nm
enpr_A=1.76;          %n del vidrio (sólo hay parte real)

%Oro, 633 nm
e1n_A=0.1726;         %Parte real de n (oro)
e1k_A=3.4218;         %Parte imaginaria de n (oro)
e1r_A=e1n_A^2-e1k_A^2;    %Parte real de epsilon (oro)
e1i_A=2*e1n_A*e1k_A;      %Parte imaginaria de epsilon (oro)
e1_A=complex(e1r_A,e1i_A); %Epsilon compleja (oro)
d1_A=47e-9;           %Grosor de la capa de oro (m)

%Agua
n2_A=1.3342;               %n del analito (sólo hay parte real)
e2_A=n2_A^2;            %epsilon del analito (sólo hay parte real)

%---ESCANEO DE LOS ÁNGULOS Y SOLUCIÓN DE LAS ECUACIONES DE FRESNEL

for x=1:(length(theta_ext))
    theta_A=theta_ext(x);
	rpr1_A=(cos(theta_A)/enpr_A - sqrt(e1_A-(enpr_A^2)*(sin(theta_A)^2))/e1_A)/...
		(cos(theta_A)/enpr_A + sqrt(e1_A-(enpr_A^2)*(sin(theta_A)^2))/e1_A);
	r12_A=(sqrt(e1_A-(enpr_A^2)*(sin(theta_A)^2))/e1_A -...
		sqrt(e2_A-(enpr_A^2)*(sin(theta_A)^2))/e2_A) /...
		(sqrt(e1_A-(enpr_A^2)*(sin(theta_A)^2))/e1_A +...
		sqrt(e2_A-(enpr_A^2)*(sin(theta_A)^2))/e2_A);
	kz1d1_A=omega/c*d1_A*sqrt(e1_A-(enpr_A^2)*(sin(theta_A)^2));
	rpr12_A=(rpr1_A+ r12_A*exp(2*i*kz1d1_A))/(1+rpr1_A*r12_A*exp(2*i*kz1d1_A));
	ref_A(x)=rpr12_A*(conj(rpr12_A))*100;
end

refmin_A=find(ref_A==min(ref_A));
thetacr_A=theta_deg(refmin_A); %Posición del ángulo crítico (mín. reflectancia)

% DEFINIMOS MATRIZ QUE RESTA LAS INTENSIDADES OBTENIDAS
resta=ref_A-ref

%---DIBUJO DE LA GRÁFICA
clf
%plot(theta_deg,resta, '--')

%Para que aparezcan juntas (mal hecho)
yyaxis left
plot(theta_deg,ref, 'k','DisplayName',...
    ['n = 1.333,       \theta_{min} = ', num2str(thetacr), '°'])
hold on;
plot(theta_deg,ref_A, 'c','DisplayName',...
    ['n = 1.3342,     \theta_{min} = ', num2str(thetacr_A), '°'])
hold off;

yyaxis right
plot(theta_deg,resta)

yyaxis left
xlabel('Ángulo \theta (°)')
ylabel('% Reflectancia')

yyaxis right
ylabel('\Delta% Reflectancia')