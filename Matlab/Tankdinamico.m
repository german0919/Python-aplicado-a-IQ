clear 
clc
%Datos
Fe=0.1;
Fs=0;
D=1;
A=(pi/4)*D^2;
%Planteo valores iniciales
h(1)=1; %Al vector h en la posicion 1 tome el valor de 1
t(1)=0; %Al vector t en la posicion 1 tome el valor de 0
t0=0;
tf=60;
n=10; %Numero de veces que voy a iterar el metodo de euler
dt=(tf-t0)/n; %Tamaño de paso
%Metodo de euler
for i=1:n
    f=(Fe-Fs)/A;
    h(i+1)=h(i)+dt*f;
    t(i+1)=t(i)+dt;
end 
plot(t,h);
xlabel('Tiempo (min)');
ylabel('Altura del líquido (m)');
title('Comportamiento dinámico de la altura del líquido dentro del tanque');
grid

    


