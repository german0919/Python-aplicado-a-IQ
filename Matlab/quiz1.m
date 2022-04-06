clear
clc
H=1.8;
Dt=2;
Do=0.05;
At=(pi/4)*Dt^2;
Ao=(pi/4)*Do^2;
h(1)=((0.01/Ao)^2)/(2*9.8);
t(1)=0;
t0=0;
tf=3600;
n=1000;
dt=(tf-t0)/n;
for i=1:n
    tm(i)=t(i)/60;
    Fe(i)=(1.67e-4)*tm(i)+0.01;
    f=(Fe(i)-Ao*sqrt(2*9.8*h(i)))/At;
    h(i+1)=h(i)+dt*f;
    t(i+1)=t(i)+dt;
end
t=t/60;
plot(t,h);
xlabel('Tiempo (min)');
ylabel('Altura de líquido (m)');
title('Comportamiento de la altura respecto al tiempo');
grid
k=1;
while h(k)<H
        k=k+1;
end
disp('El tiempo que tarda en alcanzar la altura fijada es (en min):');
disp(t(k));
