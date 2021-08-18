clear all;
figure(3)
clf
%---PARÁMETROS DE LAS DISTINTAS CAPAS
%Vidrio (SF10), 633 nm
enpr=1.76;          %n del vidrio (sólo hay parte real)

%Oro, 633 nm
e1n=0.1726;         %Parte real de n (oro)
e1k=3.4218;         %Parte imaginaria de n (oro)
e1r=e1n^2-e1k^2;    %Parte real de epsilon (oro)
e1i=2*e1n*e1k;      %Parte imaginaria de epsilon (oro)
e1=complex(e1r,e1i); %Epsilon compleja (oro)
d1=45e-9;           %Grosor de la capa de oro (m)

%Agua y etanol
n2=1.333;             %n del agua (sólo hay parte real)
n3=1.336;             %n del agua + 5% de etanol (sólo hay parte real)
n4=1.3395;            %n del agua + 10% de etanol (sólo hay parte real)
n5=1.3440;            %n del agua + 16% de etanol (sólo hay parte real)
n6=1.3469;            %n del agua + 20% de etanol (sólo hay parte real)
n7=1.3511;            %n del agua + 26% de etanol (sólo hay parte real)
n8=1.4201;            %n de agua + 50% de sacarosa (sólo hay parte real)

n9=1.3606;            %n de agua + 18% de sacarosa (sólo hay parte real)
n10=1.3706;           %n de agua + 24% de sacarosa (sólo hay parte real)
n11=1.3812;           %n de agua + 30% de sacarosa (sólo hay parte real)
n12=1.3922;           %n de agua + 36% de sacarosa (sólo hay parte real)
n13=1.4038;           %n de agua + 42% de sacarosa (sólo hay parte real)
n14=1.4118;           %n de agua + 46% de sacarosa (sólo hay parte real)

e2=n2^2;              %epsilon del agua (sólo hay parte real)
e3=n3^2;              %epsilon del agua + 5% de etanol (sólo hay parte real)
e4=n4^2;              %epsilon del agua + 10% de etanol (sólo hay parte real)
e5=n5^2;              %epsilon del agua + 16% de etanol (sólo hay parte real)
e6=n6^2;              %epsilon del agua + 20% de etanol (sólo hay parte real)
e7=n7^2;              %epsilon del agua + 26% de etanol (sólo hay parte real)
e8=n8^2;              %epsilon de agua + 50% de sacarosa (sólo hay parte real)

e9=n9^2;              %epsilon de agua + 18% de sacarosa (sólo hay parte real)
e10=n10^2;            %epsilon de agua + 24% de sacarosa (sólo hay parte real)
e11=n11^2;            %epsilon de agua + 30% de sacarosa (sólo hay parte real)
e12=n12^2;            %epsilon de agua + 36% de sacarosa (sólo hay parte real)
e13=n13^2;            %epsilon de agua + 42% de sacarosa (sólo hay parte real)
e14=n14^2;            %epsilon de agua + 46% de sacarosa (sólo hay parte real)

e=[e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14];

%---CONSTANTES Y PARÁMETROS DADOS 
c=2.99792458e8;
lambda=632.8e-9;
omega=2*pi/lambda*c;

%---ESCANEO DE LOS ÁNGULOS Y SOLUCIÓN DE LAS ECUACIONES DE FRESNEL
ang0=45;            %Límite inferior
ang1=70;            %Límite superior
vals=10000;
interval=ang1-ang0;
angmat=ang0:(interval/vals):ang1; %Vector de ángulos

%for ii=1:(length(e))
    for x=1:(length(angmat))
        theta=(ang0+((x-1)*interval)/vals)/180*pi;
        rpr1=(cos(theta)/enpr - sqrt(e1-(enpr^2)*(sin(theta)^2))/e1)/...
            (cos(theta)/enpr + sqrt(e1-(enpr^2)*(sin(theta)^2))/e1);
        r12=(sqrt(e1-(enpr^2)*(sin(theta)^2))/e1 -...
            sqrt(e2 -(enpr^2)*(sin(theta)^2))/e2) /...
            (sqrt(e1-(enpr^2)*(sin(theta)^2))/e1 +...
            sqrt(e2-(enpr^2)*(sin(theta)^2))/e2);
        kz1d1=omega/c*d1*sqrt(e1-(enpr^2)*(sin(theta)^2));
        rpr12=(rpr1+ r12*exp(2*i*kz1d1))/(1+rpr1*r12*exp(2*i*kz1d1));
        ref1(x)=rpr12*(conj(rpr12))*100;
    end
%    eval(['ref' num2str(ii) '=ref(x);']);
%end

for x=1:(length(angmat))
	theta1=(ang0+((x-1)*interval)/vals)/180*pi;
	rpr3=(cos(theta1)/enpr - sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1)/...
		(cos(theta1)/enpr + sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1);
	r34=(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 -...
		sqrt(e3-(enpr^2)*(sin(theta1)^2))/e3) /...
		(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 +...
		sqrt(e3-(enpr^2)*(sin(theta1)^2))/e3);
	kz2d2=omega/c*d1*sqrt(e1-(enpr^2)*(sin(theta1)^2));
	rpr34=(rpr3+ r34*exp(2*i*kz2d2))/(1+rpr3*r34*exp(2*i*kz2d2));
	ref2(x)=rpr34*(conj(rpr34))*100;
end

for x=1:(length(angmat))
	theta1=(ang0+((x-1)*interval)/vals)/180*pi;
	rpr3=(cos(theta1)/enpr - sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1)/...
		(cos(theta1)/enpr + sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1);
	r34=(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 -...
		sqrt(e4-(enpr^2)*(sin(theta1)^2))/e4) /...
		(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 +...
		sqrt(e4-(enpr^2)*(sin(theta1)^2))/e4);
	kz2d2=omega/c*d1*sqrt(e1-(enpr^2)*(sin(theta1)^2));
	rpr34=(rpr3+ r34*exp(2*i*kz2d2))/(1+rpr3*r34*exp(2*i*kz2d2));
	ref3(x)=rpr34*(conj(rpr34))*100;
end

for x=1:(length(angmat))
	theta1=(ang0+((x-1)*interval)/vals)/180*pi;
	rpr3=(cos(theta1)/enpr - sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1)/...
		(cos(theta1)/enpr + sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1);
	r34=(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 -...
		sqrt(e5-(enpr^2)*(sin(theta1)^2))/e5) /...
		(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 +...
		sqrt(e5-(enpr^2)*(sin(theta1)^2))/e5);
	kz2d2=omega/c*d1*sqrt(e1-(enpr^2)*(sin(theta1)^2));
	rpr34=(rpr3+ r34*exp(2*i*kz2d2))/(1+rpr3*r34*exp(2*i*kz2d2));
	ref4(x)=rpr34*(conj(rpr34))*100;
end

for x=1:(length(angmat))
	theta1=(ang0+((x-1)*interval)/vals)/180*pi;
	rpr3=(cos(theta1)/enpr - sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1)/...
		(cos(theta1)/enpr + sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1);
	r34=(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 -...
		sqrt(e6-(enpr^2)*(sin(theta1)^2))/e6) /...
		(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 +...
		sqrt(e6-(enpr^2)*(sin(theta1)^2))/e6);
	kz2d2=omega/c*d1*sqrt(e1-(enpr^2)*(sin(theta1)^2));
	rpr34=(rpr3+ r34*exp(2*i*kz2d2))/(1+rpr3*r34*exp(2*i*kz2d2));
	ref5(x)=rpr34*(conj(rpr34))*100;
end

for x=1:(length(angmat))
	theta1=(ang0+((x-1)*interval)/vals)/180*pi;
	rpr3=(cos(theta1)/enpr - sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1)/...
		(cos(theta1)/enpr + sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1);
	r34=(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 -...
		sqrt(e7-(enpr^2)*(sin(theta1)^2))/e7) /...
		(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 +...
		sqrt(e7-(enpr^2)*(sin(theta1)^2))/e7);
	kz2d2=omega/c*d1*sqrt(e1-(enpr^2)*(sin(theta1)^2));
	rpr34=(rpr3+ r34*exp(2*i*kz2d2))/(1+rpr3*r34*exp(2*i*kz2d2));
	ref6(x)=rpr34*(conj(rpr34))*100;
end

for x=1:(length(angmat))
	theta1=(ang0+((x-1)*interval)/vals)/180*pi;
	rpr3=(cos(theta1)/enpr - sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1)/...
		(cos(theta1)/enpr + sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1);
	r34=(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 -...
		sqrt(e8-(enpr^2)*(sin(theta1)^2))/e8) /...
		(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 +...
		sqrt(e8-(enpr^2)*(sin(theta1)^2))/e8);
	kz2d2=omega/c*d1*sqrt(e1-(enpr^2)*(sin(theta1)^2));
	rpr34=(rpr3+ r34*exp(2*i*kz2d2))/(1+rpr3*r34*exp(2*i*kz2d2));
	ref7(x)=rpr34*(conj(rpr34))*100;
end

for x=1:(length(angmat))
	theta1=(ang0+((x-1)*interval)/vals)/180*pi;
	rpr3=(cos(theta1)/enpr - sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1)/...
		(cos(theta1)/enpr + sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1);
	r34=(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 -...
		sqrt(e9-(enpr^2)*(sin(theta1)^2))/e9) /...
		(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 +...
		sqrt(e9-(enpr^2)*(sin(theta1)^2))/e9);
	kz2d2=omega/c*d1*sqrt(e1-(enpr^2)*(sin(theta1)^2));
	rpr34=(rpr3+ r34*exp(2*i*kz2d2))/(1+rpr3*r34*exp(2*i*kz2d2));
	ref8(x)=rpr34*(conj(rpr34))*100;
end

for x=1:(length(angmat))
	theta1=(ang0+((x-1)*interval)/vals)/180*pi;
	rpr3=(cos(theta1)/enpr - sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1)/...
		(cos(theta1)/enpr + sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1);
	r34=(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 -...
		sqrt(e10-(enpr^2)*(sin(theta1)^2))/e10) /...
		(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 +...
		sqrt(e10-(enpr^2)*(sin(theta1)^2))/e10);
	kz2d2=omega/c*d1*sqrt(e1-(enpr^2)*(sin(theta1)^2));
	rpr34=(rpr3+ r34*exp(2*i*kz2d2))/(1+rpr3*r34*exp(2*i*kz2d2));
	ref9(x)=rpr34*(conj(rpr34))*100;
end

for x=1:(length(angmat))
	theta1=(ang0+((x-1)*interval)/vals)/180*pi;
	rpr3=(cos(theta1)/enpr - sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1)/...
		(cos(theta1)/enpr + sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1);
	r34=(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 -...
		sqrt(e11-(enpr^2)*(sin(theta1)^2))/e11) /...
		(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 +...
		sqrt(e11-(enpr^2)*(sin(theta1)^2))/e11);
	kz2d2=omega/c*d1*sqrt(e1-(enpr^2)*(sin(theta1)^2));
	rpr34=(rpr3+ r34*exp(2*i*kz2d2))/(1+rpr3*r34*exp(2*i*kz2d2));
	ref10(x)=rpr34*(conj(rpr34))*100;
end

for x=1:(length(angmat))
	theta1=(ang0+((x-1)*interval)/vals)/180*pi;
	rpr3=(cos(theta1)/enpr - sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1)/...
		(cos(theta1)/enpr + sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1);
	r34=(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 -...
		sqrt(e12-(enpr^2)*(sin(theta1)^2))/e12) /...
		(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 +...
		sqrt(e12-(enpr^2)*(sin(theta1)^2))/e12);
	kz2d2=omega/c*d1*sqrt(e1-(enpr^2)*(sin(theta1)^2));
	rpr34=(rpr3+ r34*exp(2*i*kz2d2))/(1+rpr3*r34*exp(2*i*kz2d2));
	ref11(x)=rpr34*(conj(rpr34))*100;
end

for x=1:(length(angmat))
	theta1=(ang0+((x-1)*interval)/vals)/180*pi;
	rpr3=(cos(theta1)/enpr - sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1)/...
		(cos(theta1)/enpr + sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1);
	r34=(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 -...
		sqrt(e13-(enpr^2)*(sin(theta1)^2))/e13) /...
		(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 +...
		sqrt(e13-(enpr^2)*(sin(theta1)^2))/e13);
	kz2d2=omega/c*d1*sqrt(e1-(enpr^2)*(sin(theta1)^2));
	rpr34=(rpr3+ r34*exp(2*i*kz2d2))/(1+rpr3*r34*exp(2*i*kz2d2));
	ref12(x)=rpr34*(conj(rpr34))*100;
end

for x=1:(length(angmat))
	theta1=(ang0+((x-1)*interval)/vals)/180*pi;
	rpr3=(cos(theta1)/enpr - sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1)/...
		(cos(theta1)/enpr + sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1);
	r34=(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 -...
		sqrt(e14-(enpr^2)*(sin(theta1)^2))/e14) /...
		(sqrt(e1-(enpr^2)*(sin(theta1)^2))/e1 +...
		sqrt(e14-(enpr^2)*(sin(theta1)^2))/e14);
	kz2d2=omega/c*d1*sqrt(e1-(enpr^2)*(sin(theta1)^2));
	rpr34=(rpr3+ r34*exp(2*i*kz2d2))/(1+rpr3*r34*exp(2*i*kz2d2));
	ref13(x)=rpr34*(conj(rpr34))*100;
end

refmin1=find(ref1==min(ref1));
thetacr1=angmat(refmin1); %Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico 1 es: ');
disp(thetacr1);

refmin2=find(ref2==min(ref2));
thetacr2=angmat(refmin2); %Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico 2 es: ');
disp(thetacr2);

refmin3=find(ref3==min(ref3));
thetacr3=angmat(refmin3); %Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico 3 es: ');
disp(thetacr3);

refmin4=find(ref4==min(ref4));
thetacr4=angmat(refmin4); %Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico 4 es: ');
disp(thetacr4);

refmin5=find(ref5==min(ref5));
thetacr5=angmat(refmin5); %Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico 5 es: ');
disp(thetacr5);

refmin6=find(ref6==min(ref6));
thetacr6=angmat(refmin6); %Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico 6 es: ');
disp(thetacr6);

refmin7=find(ref7==min(ref7));
thetacr7=angmat(refmin7); %Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico 7 es: ');
disp(thetacr7);

refmin8=find(ref8==min(ref8));
thetacr8=angmat(refmin8); %Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico 8 es: ');
disp(thetacr8);

refmin9=find(ref9==min(ref9));
thetacr9=angmat(refmin9); %Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico 9 es: ');
disp(thetacr9);

refmin10=find(ref10==min(ref10));
thetacr10=angmat(refmin10); %Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico 10 es: ');
disp(thetacr10);

refmin11=find(ref11==min(ref11));
thetacr11=angmat(refmin11); %Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico 11 es: ');
disp(thetacr11);

refmin12=find(ref12==min(ref12));
thetacr12=angmat(refmin12); %Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico 12 es: ');
disp(thetacr12);

refmin13=find(ref13==min(ref13));
thetacr13=angmat(refmin13); %Posición del ángulo crítico (mín. reflectancia)
disp('El angulo critico 13 es: ');
disp(thetacr13);

%---DIBUJO DE LA GRÁFICA
plot(angmat,ref1,'DisplayName','n = 1.333');
hold on;
plot(angmat,ref2,'DisplayName','n = 1.336');
plot(angmat,ref3,'DisplayName','n = 1.3395');
plot(angmat,ref4,'DisplayName','n = 1.3440');
plot(angmat,ref5,'DisplayName','n = 1.3469');
plot(angmat,ref6,'DisplayName','n = 1.3511');
plot(angmat,ref7,'DisplayName','n = 1.3606');
plot(angmat,ref8,'DisplayName','n = 1.3706');
plot(angmat,ref9,'DisplayName','n = 1.3812');
plot(angmat,ref10,'DisplayName','n = 1.3922');
plot(angmat,ref11,'DisplayName','n = 1.4038');
plot(angmat,ref12,'DisplayName','n = 1.4118');
plot(angmat,ref13,'DisplayName','n = 1.4201');
hold off;

title('SIMULACIÓN DE CURVA SPR A DIFERENTES RI')
ylabel('Reflectancia (%)')
xlabel('Ángulo \theta (°)')