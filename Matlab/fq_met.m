clc
clear all

Cao=10; %Concentración inicial de A [M]
Cbo=0; %Concentración inicial de B [M]
Cco=0; %Concentración inicial de C [M]

%Constantes cineticas de las reacciones reversibles 

K1=input('Ingrese el valor de la constante cinetica de la reaccion directa 1 [h-1] ') ; %Constante cinetica directa 1 (A->B) [h-1]
K11=input('Ingrese el valor de la constante cinetica de la reaccion inversa 1 [h-1] ') ; %Constante cinetica inversa -1 (B->A) [h-1]
K2=input('Ingrese el valor de la constante cinetica de la reaccion directa 2 [h-1] ') ; %Constante cinetica directa 2 [h-1]
K22=input('Ingrese el valor de la constante cinetica de la reaccion inversa 2 [h-1] ') ; %Constante cinetica inversa -2 [h-1]
Xa=input('Ingrese la conversi�n deseada de A ');

        n=1;
        X(n)=0; t(n)=0;Ca(n)=Cao; Cb(n)=Cbo; Cc(n)=Cco; deltat= 0.01;
        Cai(n)=0; Cbi(n)=0; Cci(n)=0;
        while X(n)< Xa
                    
        %Ecuaciones diferenciales
       
        dCadt=-K1*Ca(n)+K11*Cb(n);
        dCbdt=K1*Ca(n)+K22*Cc(n)-Cb(n)*(K11+K2);
        dCcdt=K2*Cb(n)-K22*Cc(n);
        
       
        %Metodo de Euler modificado
        Cai(n)=Ca(n)+dCadt*(deltat)*0.5;
        Cbi(n)=Cb(n)+dCbdt*(deltat)*0.5;
        Cci(n)=Cc(n)+dCcdt*(deltat)*0.5;
        
        %Derivada evaluada en el punto medio
        dCadti=-K1*Cai(n)+K11*Cbi(n);
        dCbdti=K1*Cai(n)+K22*Cci(n)-Cbi(n)*(K11+K2);
        dCcdti=K2*Cbi(n)-K22*Cci(n);
                
        
        Ca(n+1)=Ca(n)+dCadti*deltat;
        Cb(n+1)=Cb(n)+dCbdti*deltat;
        Cc(n+1)=Cc(n)+dCcdti*deltat;
        
                
        t(n+1)=t(n)+deltat;
        X(n+1)=(Cao-Ca(n+1))/Cao;
       
        n=n+1;
        end
        %Gráficas
        
        subplot(2,2,1)
        plot(t,Ca,t,Cb,t,Cc)
        xlabel('Tiempo (h)')
        title('Concentración vs tiempo')
        legend('Concentración A','Concentración B','Concentración C')
        
        subplot(2,2,2)
        plot(t,Ca,'b')
        xlabel('Tiempo (h)')
        ylabel('Concentración A')
        title('Concentración A vs tiempo')
        
        
        subplot(2,2,3)
        plot(t,Cb,'r')
        xlabel('Tiempo (h)')
        ylabel('Concentración B')
        title('Concentración B vs tiempo')
      
        
        subplot(2,2,4)
        plot(t,Cc,'y')
        xlabel('Tiempo (h)')
        ylabel('Concentración C')
        title('Concentración C vs tiempo')
       
        
        
  

 










