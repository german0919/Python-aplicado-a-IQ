clear 
clc
%German Hernández-2171842
%Datos
%---------------------------------
rho=998; %kg/m^3
W=8000; %W
nb=0.70;
mu=1.002e-3; %kg/m
e=0.000045;
g=9.81;
L1=36;%m
L2=36; %m
D1=0.04;%m
D2=0.08; %m
za=5; %m
zb=13; %m
A1=(pi/4)*(D1^2); %m^2
A2=(pi/4)*(D2^2);

%---------------------------------
er11=10;
er22=10;
i=1;
%Supongo v
v1(i)=5;
v2(i)=10;
while er11>1e-4 &&  er22>1e-4
    Q1=v1(i)*A1;
    Q2=v2(i)*A2;
    Re1=(rho*v1(i)*D1)/mu;
    Re2=(rho*v2(i)*D2)/mu;
    Qt=Q1+Q2;
    
    %Calculo f1
    k=1;
    x1(k)=0.01;
    er12=10;
    while er12>0.01
        x1(k+1)=1.14-2*log10((e/D1)+(9.35/Re1)*x1(k));
        er12=(abs(x1(k+1)-x1(k))/x1(k+1))*100;
        k=k+1;
    end
    ff1=(1/x1(k))^2;
    %Calculo f2
    m=1;
    x2(m)=0.01;
    er21=10;
    while er21>0.01
        x2(m+1)=1.14-2*log10((e/D2)+(9.35/Re2)*x2(m));
        er21=(abs(x2(m+1)-x2(m))/x2(m+1))*100;
        m=m+1;
    end
    ff2=(1/x2(m))^2;
    %Calculo hb
    hb=(W*nb)/(Qt*g*rho);
    %Calculo hl/recordando que hL1=hL2
    hl=hb-(zb-za);
    %Compruebo v
    v1(i+1)=((hl*D1*2*g)/(ff1*L1))^0.5;
    v2(i+1)=((hl*D2*2*g)/(ff2*L2))^0.5;
    er11=(abs(v1(i+1)-v1(i))/v1(i+1))*100;
    er22=(abs(v2(i+1)-v2(i))/v2(i+1))*100;
    i=i+1;
end
disp('Las velocidades son:');disp(v1(i)); disp(v2(i))
disp('El caudal total es:');disp(Qt);

