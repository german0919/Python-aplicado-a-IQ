clear 
clc
m1ac=0.2*1000;
m1ag=0.8*1000;
m4ag=m1ag;
m4b=0.07*m4ag;
%Debemos calcular 0.07*0.8*1000 porque m2 debe ser mayor que esa cantidad
%Ahora si se supone el valor
%Vamos a colocar un vector para hacer el problema
m2=100:1:3000; %valores de 100 a 3000 de 1 en 1
n=length(m2);
for i=1:n
m3b(i)=m2(i)-m4b;
m3ac(i)=m1ac/(1+((m4ag+m4b)/(4*m3b(i))));
U(i)=3*m3ac(i)-0.03*m2(i);
end
plot(m2,U);
xlabel('Flujo de benceno [kg/h]')
ylabel('Utilidad [$]')
title('Utilidad vs flujo de agente extractor')
grid
Um=max(U);
k=1;
while Um>U(k)
    k=k+1;
end
disp('El flujo óptimo de benceno es (en kg/s):')
disp(m2(k));
disp('La utilidad máxima es (en $):')
disp(U(k))
    

