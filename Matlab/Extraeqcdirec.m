clear
clc
ye=0.05;
ys=0.005;
xe=0;
ER=1.1;
m=0.71;
i=1;
y(i)=ye;
x(i)=xe;
%mientras que la composicion de salida sea mayor que la fraccion  en el
%refinado haga el calculo
while y(i)>ys
    x(i+1)=(y(i)+x(i)*ER/(m+ER));
    y(i+1)=m*x(i+1);
    i=i+1;
    if i>100
        disp('Etapas infinitas. El flujo de agente extractor no es suficiente para la separación deseada');
        break
    end    
end
if i<50
   disp('El número de etapas es :');
   disp(i-1);
end 