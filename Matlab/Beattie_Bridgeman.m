clear all
clc
close all
%P=kpa T=K V=m^3
P=405.3; 
T=-10+273.15; 
%Parametros para el Co2
Ru=8.314;%kPa/m^3*kmol
A0=507.2836;
a=0.07132;
B0=0.10476;
b=0.07235;
c=6.60e5;
%volumen inicial
v=Ru*T/P;
%Densidad inicial
rho=1/v;
%Para hallar la derivada en la ventana de comandos
%syms Ru T B0 b A0 a P c rho
%(Ru*T*rho^2)(1-((rho*c)/T^3))(1/rho + B0 - B0*b*rho)- (A0 - A0*a)*rho^3 - P
i=1;
tol=1e-3;
er=tol+1;
while er<tol
    f=(Ru*T*rho^2)*(1-((rho*c)/T^3))*(1/rho+B0-B0*b*rho)-(A0-A0*a)*rho^3-P;
    df =Ru*T*rho^2*(B0*b + 1/rho^2)*((c*rho)/T^3 - 1) - 3*rho^2*(A0 - A0*a) - (Ru*c*rho^2*(B0 + 1/rho - B0*b*rho))/T^2 - 2*Ru*T*rho*((c*rho)/T^3 - 1)*(B0 + 1/rho - B0*b*rho); 
    rhon= rho-f/df;
    er=abs((rhon-rho)/rhon);
    fprintf('%d\t %f\t %f\f %f\n',i,rho,rhon,er)
    rho=rhon;
    i=i+1;
end

disp('El valor de la densidad es :')
disp(rho)


