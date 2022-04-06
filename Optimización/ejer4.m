clear
clc
func= @(x) -(3*x(1)*exp(-0.1*x(1)*x(6))+4*x(2)+x(3).^2+7*x(4)+(10/x(5))+x(6));
Aeq=[0 0 0 -1 1 -1;1 1 1 1 1 1] ; beq=[0.1;10];
lb=[-20;-20;-20;-20;-20;-20]; ub=[20;20;20;20;20;20]; 
A=[-2 -1 -1 -3 0 0;8 3 4 -1 1 0;2 6 1 3 0 1;1 4 5 2 0 0]; b=[-2;10;13;18];
nonlcon=[];options=[];
x0=[10 10 12 17 10 15];
[x,fval,exitflag]=fmincon(func,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
x
