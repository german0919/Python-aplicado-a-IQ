clear 
clc
ye=0.05; %comp entrada del refinado
ys=0.005; %comp de salida del soluto
xe=0; %fraccion de entrada del extracto
m=0.71;
ER=0.72:0.01:1.5;%relacion extracto refinado
n=length(ER);
for k=1:n
    i=1; %contador
    y(i)=ys; %valor inicial del vector y
    x(i)=xe; %valor inciial del vector x
    %Calculo etapa por etapa
    while y(i)<ye
        x(i+1)=y(i)/m;
        y(i+1)=y(i)-ER(k)*(x(i+1)-x(i));
        i=i+1;
    end
    et(k)=i-1; %guardar el numero de etapas para cada posicion en ER(k )
end
plot(ER,et,'linewidth',1.5)
xlabel('Flujo de extracto por mol de gas libre de soluto');
ylabel('Número de etapas ideales');
title('Etapas vs. flujo de agente extractor','fontsize',15);
grid
axis([0.71 1.5, 0 et(1)])