#Ejercicio 2- GUIA 2D
clc;
clear all;
xnode=[0 0;
       0.5 0;
       1 0;
       2.5 0;
       0 1;
       0.5 1;
       1 1];
       
icone=[1 2 6 5;
       2 3 7 6;
       3 4 7 -1];
model.nnodes=7;
model.nelem=3;
model.kx=1;
model.ky=1;
model.c=0;
G=@(x,y) 10.*x;
model.G=G(xnode(:,1),xnode(:,2));
model.ts=-1;
model.rho=0;
model.cp=0;
model.maxit=5000;
model.tol=0.00000000001;
model.dt=0.01;
model.PHI_n=[0 0 0 0 0 0 0 0 0]';

DIR=[1 0;2 0;3 0;4 0;5 0];
NEU=[5 6 0;6 7 0;4 7 0];
ROB=[];
PUN=[];

[PHI,Q]= fem2d_heat(xnode,icone,DIR,NEU,ROB,PUN,model);
fem2d_heat_graph_mesh(PHI,Q,xnode,icone,0,1);