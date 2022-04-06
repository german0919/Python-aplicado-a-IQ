clear
clc
func= @(x) (x(1)-3).^2+(x(2)+4).^2+exp(5*x(3));
Aeq=[] ; beq=[];
lb=zeros(1,3); ub=[]; 
A=[1 1 1]; b=1;
nonlcon=[];options=[];
x0=[1 1 1];
[x,fval,exitflag,output,lambda,grad,hessian]=fmincon(func,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
x