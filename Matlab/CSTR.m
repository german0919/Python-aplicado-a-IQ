%German Hernández-2171842
clear 
clc
%Datos
F=100;
ce=12*1.2; %debia multiplicar por 1.2 dado que me preguntaban que sucedia si cambiaba la variable de entrada
k=52; %h
V=5;
%Condiciones iniciales
cs(1)=(100*12)/(100+k*V); %estado estable
t(1)=0;
t0=0;
tf=0.2;
n=1e6;
dt=(tf-t0)/n; %Tamaño de paso
%Metodo de euler
for i=1:n
    fu(i)=(F*ce-F*cs(i)-k*cs(i)*V)/V;
    cs(i+1)=cs(i)+dt*fu(i);
    t(i+1)=t(i)+dt;
end 
plot(t,cs);
xlabel('Tiempo (h)');
ylabel('Cs (Kmol/m^3)');
title('CSTR');
grid

