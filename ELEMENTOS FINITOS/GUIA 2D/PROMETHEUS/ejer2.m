#Ejercicio 2- GUIA 2D
clc;
clear all;
xnode=[0 0;
       0.5 0;
       1 0;
       0 0.5;
       0.5 0.5;
       1 0.5;
       0 1;
       0.5 1;
       1 1];
icone=[1 2 5 4;
       2 3 6 5;
       4 5 8 7;
       5 6 9 8];

model.nnodes=9;
model.nelem=4;
model.kx=1;
model.ky=1;
model.c=0;
G=@(x,y) 100.*((x-0.5).^2 + (y-0.5).^2);
model.G=G(xnode(:,1),xnode(:,2));
model.ts=2;
model.rho=1;
model.cp=1;
model.maxit=1000;
model.tol=0.0001;
model.dt=1e-3;
model.PHI_n=[0 0 0 0];

DIR=[1 0;
     2 0;
     3 0;
     4 0;
     5 0;
     6 0;
     7 0;
     8 0;
     9 0;];
NEU=[];
ROB=[];
PUN=[];

[PHI,Q]= fem2d_heat(xnode,icone,DIR,NEU,ROB,PUN,model);
fem2d_heat_graph_mesh(PHI,Q,xnode,icone,0,0);
