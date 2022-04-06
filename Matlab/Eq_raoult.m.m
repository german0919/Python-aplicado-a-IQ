clear
clc
disp('GRÁFICA DE EQUILIBRIO LÍQUIDO/VAPOR PARA UN SISTEMA BINARIO MODELADO POR LA LEY DE RAOULT')
P=input('Ingrese el valor de la presión en kPa: ');
T1=Antoine12(P);
T2=Antoine22(P);
x=linspace(T1,T2);
n=length(x);
for k=1:n
    p1(k)=Antoine11(x(k));
    p2(k)=Antoine21(x(k));
    X(k)=(P-p2(k))/(p1(k)-p2(k));
    Y(k)=X(k)*p1(k)/P;
end
plot(X,x,'b-',Y,x,'r-');
axis([0 1,T1-5 T2+5])
xlabel('x,y');
ylabel('Temperatura');
grid
t=title('Equilibrio líquido/vapor acetonitrilo/nitrometano');
set(t,'fontsize',20);
function [y]=Antoine11(x)
A=14.2724;
B=2945.47;
C=224;
y=exp(A-B/(x+C));
end
function [y]=Antoine12(x)
A=14.2724;
B=2945.47;
C=224;
y=B/(A-log(x))-C;
end
function [y]=Antoine21(x)
A=14.2043;
B=2972.64;
C=209;
y=exp(A-B/(x+C));
end
function [y]=Antoine22(x)
A=14.2043;
B=2972.64;
C=209;
y=B/(A-log(x))-C;
end