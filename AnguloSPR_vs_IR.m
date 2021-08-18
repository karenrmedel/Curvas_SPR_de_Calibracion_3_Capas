clear all;
figure(2)
clf

%---PARÁMETROS DE LAS DISTINTAS CAPAS
%Agua y etanol
e2=1.333;             %n del agua (sólo hay parte real)
e3=1.334;             %n del agua + 5% de etanol (sólo hay parte real)
e4=1.335;            %n del agua + 10% de etanol (sólo hay parte real)
e5=1.336;            %n del agua + 16% de etanol (sólo hay parte real)
e6=1.337;            %n del agua + 20% de etanol (sólo hay parte real)
e7=1.338;
e8=1.3395;
e9=1.34;
e10=1.341;
e11=1.342;
e12=1.343;
e13=1.344;
e14=1.345;
e15=1.346;
e16=1.3469;
e17=1.348;
e18=1.349;
e19=1.350;
e20=1.3511;
e21=1.3606;
e22=1.3706;
e23=1.3812;
e24=1.3922;
e25=1.4038;
e26=1.4118;
e27=1.4201;

e=[e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20,e21,e22,e23,e24,e25,e26,e27];

%{
e=[e2,e3,e4,e5,e6,e7,e8,e9];
%}

%Ángulos críticos
thetacr2=55.46;            %angulocr del agua (sólo hay parte real)
thetacr3=55.5325;            %angulocr del agua + 5% de etanol (sólo hay parte real)
thetacr4=55.605;            %angulocr del agua + 10% de etanol (sólo hay parte real)
thetacr5=55.6775;            %angulocr del agua + 16% de etanol (sólo hay parte real)
thetacr6=55.75;            %angulocr del agua + 20% de etanol (sólo hay parte real)
thetacr7=55.8225;
thetacr8=55.93;
thetacr9=55.9675;
thetacr10=56.04;
thetacr11=56.1125;
thetacr12=56.185;
thetacr13=56.126;
thetacr14=56.3325;
thetacr15=56.4075;
thetacr16=56.4725;
thetacr17=56.555;
thetacr18=56.6275;
thetacr19=56.7025;
thetacr20=56.7850;            %angulocr de agua + 26% de sacarosa (sólo hay parte real)
thetacr21=57.5;
thetacr22=58.27;
thetacr23=59.105;
thetacr24=59.995;
thetacr25=60.9625;
thetacr26=61.6475;
thetacr27=62.375;


thetacr=[thetacr2,thetacr3,thetacr4,thetacr5,thetacr6,thetacr7,thetacr8,thetacr9,thetacr10,thetacr11,thetacr12,thetacr13,thetacr14...
    thetacr15,thetacr16,thetacr17,thetacr18,thetacr19,thetacr20,thetacr21,thetacr22,thetacr23,thetacr24,thetacr25,thetacr26,thetacr27];


%Quitar los comentarios del siguiente para crear el "closer"
%{
thetacr=[thetacr2,thetacr3,thetacr4,thetacr5,thetacr6,thetacr7,thetacr8,thetacr9];
%}

% GRAFICA DE ANGULO CRITICO VS IR
plot(e,thetacr,'-o');
title('SIMULACION: ÁNGULO SRP VS RI')
ylabel('Ángulo SPR (°)')
xlabel('Índice de Refracción (RIU)')