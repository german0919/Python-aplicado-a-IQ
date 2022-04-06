clear all
clc
close all
%Método de Runge kutta
func=@(x,y) y+2*x-x^2;
x=linspace(0,1.5,30);
y=zeros(size(x));
h=(1.5-0)/30;
y(1)=1;
for n=1:length(x)-1
    k1=func(x(n),y(n));
	k2=func(x(n)+h/2,y(n)+h*(k1/2));
	k3=func(x(n)+h/2,y(n)+h*(k2/2));
	k4=func(x(n)+h,y(n)+h*k3);
	y(n+1)=y(n)+(h*(k1+2*k2+2*k3+k4)/6);
end    
plot(x,y,'r-')
xlabel('x')
ylabel('y')
grid()
title('Método de Runge kutta')