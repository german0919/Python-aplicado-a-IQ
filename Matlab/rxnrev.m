clear
clc
%Reacción química reversible
F=100;
Ce=12*1.2;
k=52; %h
V=5;
Cs(1)=(100*12)/(100+k*V); %condicion inicial sup estado estable
t(1)=0;
t0=0;
tf=0.2;
n=1000;
dt=(tf-t0)/n;
for i=1:n
    f=(F*Ce-F*Cs(i)-k*V*Cs(i))/V;
    Cs(i+1)=Cs(i)+dt*f;
    t(i+1)=t(i)+dt;
end
p=plot(t,Cs);
xlabel('Tiempo (h)');
ylabel('Concentración de reactivo (kmol/m3)');
title('Comportamiento de la concentración de A respecto al tiempo');
grid
axis([0 tf, 3.2 4.1]); %establecer los limites para los ejes [xmin xmax ymin ymax]
set(p,'linewidth',1.5);