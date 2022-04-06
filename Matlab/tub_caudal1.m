clear
clc
d=0.93*1000;
vis=4e-3;
D=3*2.54/100;
e=4.6e-5;
er1=10;
i=1;
v(i)=1;
while er1>1
    %calculo para f1
    Re=(d*v(i)*D)/vis;
    k=1;
    x(k)=0.1;
    er2=10;
    while er2>0.1
        x(k+1)=1.14-2*log10((e/D)+(9.35/Re)*x(k));
        er2=(abs(x(k+1)-x(k))/x(k+1))*100;
        k=k+1;
    end
    ff=(1/x(k))^2;
    v(i+1)=(294/(577*ff+2.49))^0.5;
    er1=(abs(v(i+1)-v(i))/v(i+1))*100;
    i=i+1;
end
A=(pi/4)*(D^2);
Q=v(i)*A;
disp('El caudal de fluido es (en m3/s): ');
disp(Q);