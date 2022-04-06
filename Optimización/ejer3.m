clear
clc
func= @(x) (x(1)-1).^2+(x(1)-x(2)).^2+(x(2)-x(3)).^2+(x(3)-x(4)).^4+(x(4)-x(5)).^4;
Aeq=[] ; beq=[];
lb=[]; ub=[]; 
A=[]; b=[];
nonlcon=[@rest];options=[];
x0=[1 1 1 1 1];
[x,fval,exitflag]=fmincon(func,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
x

function [c,ceq]=rest(x)
c=[];
ceq1=x(1)+x(2).^2+x(3).^3-2-3*sqrt(2);
ceq2=x(2)-x(3).^2+x(4)+2-2*sqrt(2);
ceq3=x(1)*x(5)-2;
ceq=[ceq1;ceq2;ceq3];

end