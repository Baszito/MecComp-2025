#PRUEBA

clc;
clear all;
##
###Ejercicio A
xnode=linspace(0,1,4);
model.rho=0;
model.cp=0;
model.k=2;
model.c=0;
model.G=@(x) zeros(1,length(x))+100;
cb=[1 10 0;1 50 0];

elemFinitos(xnode,model,cb,[0 0 0 0],2)

####Ejercicio B
##xnode=linspace(0,2,4);
##model.p=0;
##model.cp=0;
##model.k=1;
##model.c=1;
##G=@(x) zeros(1,length(x));
##model.G=G(xnode);
##cb=[1 100 0;2 0 0];
##
##elemFinitos(xnode,model,cb,0)

#####Ejercicio C
##clc
##clear all
##xnode=linspace(1,5,4);
##model.p=0;
##model.cp=0;
##model.k=1;
##model.c=0;
####OBTENER VALORES APARTE CON cuad_gauss_c
##G=[167.901 88.888;9.8765 9.8765;88.888 167.9012];
##model.G=G;
##cb=[2 2 0;1 0 0];
##
##
##elemFinitos(xnode,model,cb,0)

####Ejercicio D
##xnode=linspace(0,1,4);
##model.p=0;
##model.cp=0;
##model.k=1;
##model.c=1;
##G=@(x) zeros(1,length(x))+50;
##model.G=G(xnode);
##cb=[1 10 0;3 0.2 50];
##
##
##elemFinitos(xnode,model,cb,0)

####Ejercicio E
##xnode=linspace(5,10,4);
##model.p=1;
##model.cp=1;
##model.k=2;
##model.c=0;
##model.G=[145.10 192.47;316.79 395.02;588.78 705.59];
##cb=[3 2 100;1 50 0];
##
##T=elemFinitos(xnode,model,cb,1)
##
##T_an=@(z) ((-z.^5)/40) + (1225.*z)/3 - 4600/3;
##vec_y=T_an(xnode);
##
##figure;
##plot(xnode,T,'b');
##hold on;
##grid on;
##plot(xnode,vec_y,'r');

##Ejercicio F
##xnode=linspace(0,1,4);
##model.p=2;
##model.cp=1;
##model.k=2;
##model.c=2;
##G=@(x) zeros(1,length(x))+75;
##model.G=G(xnode);
##cb=[1 0 0;3 2 10];
##
##T=elemFinitos(xnode,model,cb,1)
##
##T_an=@(z) (-5/4).*(e.^-(z+1)).*((e.^z)-1).*(11.*(e.^z) + 11-30*e);
##vec_y=T_an(xnode);
##
##figure;
##plot(xnode,T,'b');
##hold on;
##grid on;
##plot(xnode,vec_y,'r');

####Ejercicio G
##xnode=linspace(0,1,4);
##model.rho=1;
##model.cp=1;
##model.k=2;
##model.c=-2;
##model.G=@(x) zeros(1,length(x));
##cb=[1 50 0;2 5 0];
##
##
##T=elemFinitos(xnode,model,cb,[1 10000 0.0001 0.01])
##
##T_an=@(z) 73.2433.*sin(z)+50.*cos(z);
##vec_y=T_an(xnode);
##
##figure;
##plot(xnode,T,'b');
##hold on;
##grid on;
##plot(xnode,vec_y,'r');
