%Proceso de calentamiento
clear
clc
Fe=0.5;
Fv=12*1.2;
Te=90;
V=30;
d=870;
dv=1/1.694;
Cp=2.09;
L=2257;
Ts(1)=Te+(12*dv*L)/(0.5*d*Cp);
t(1)=0;
t0=0;
tf=10*60;
n=100;
dt=(tf-t0)/n;
for i=1:n
    f=(Fe*d*Cp*(Te-Ts(i))+(Fv*dv*L))/(d*V*Cp);
    Ts(i+1)=Ts(i)+dt*f;
    t(i+1)=t(i)+dt;
end
t=t/60;
tf=tf/60;
p=plot(t,Ts);
xlabel('Tiempo (min)');
ylabel('Temperatura del crudo (°C)');
title('Comportamiento de la temperatura respecto al tiempo');
grid
axis([0 tf, 106 112]);
set(p,'linewidth',1.5);


