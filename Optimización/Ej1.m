clear
clc
func= @(x) (x(1)-x(2)).^2+(x(2)-x(3)).^2+(x(3)-x(4)).^2+(x(4)-x(5)).^2;
Aeq=[1 2 3 0 0;0 1 2 3 0;0 0 1 2 3 ] ; beq=[6;6;6];
lb=[]; ub=[]; 
A=[]; b=[];
nonlcon=[];options=[];
x0=[1 2 1 3 1];
[x,fval,exitflag]=fmincon(func,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
x