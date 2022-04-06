clear 
clc
ye=0.05; %comp entrada del refinado
ys=0.005; %comp de salida del soluto
xe=0; %fraccion de entrada del extracto
m=0.71;
ER=1.1;%relacion extracto refinado
i=1; %contador
y(i)=ys; %valor inicial del vector y
x(i)=xe; %valor inciial del vector x
%Calculo etapa por etapa
while y(i)<ye
    x(i+1)=y(i)/m;
    y(i+1)=y(i)-ER*(x(i+1)-x(i));
    i=i+1;
end
disp('El número de etapas requeridas es :');
disp(i-1);