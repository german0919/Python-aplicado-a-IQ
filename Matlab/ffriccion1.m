clear 
clc
%------------------------------
%Cálculo del factor de fricción
%------------------------------
D=0.1;
e=5e-4;
Re=5e6;
i=1; %contador
er=10;
%1/sqrt(f)=x
x(i)=1; %valor supuesto(cualquiera)/con cambio de variable para disminuir el error
%Comienza el método iterativo
while er>0.1
    x(i+1)=1.14-2*log10((e/D)+(9.35/Re)*x(i)); % con la posicion i del vector x calculo x de la posicion i+1 hasta alcanzar convergencia
    er=(abs(x(i+1)-x(i))/x(i+1))*100; %dif del recalculado y el supuesto dividido en el recalculado
    i=i+1;
end
%Despejando x
ff=(1/x(i))^2; %se pone i porque es el último valor que toma i en el calculo
disp('El factor de fricción es: ');
disp(ff);

