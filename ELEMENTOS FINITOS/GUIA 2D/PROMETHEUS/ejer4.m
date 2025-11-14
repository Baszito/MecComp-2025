#Ejercicio 2- GUIA 2D
clc;
clear all;
xnode=[0 0;
       1 0
       2 0;
       0 0.5;
       1 0.5;
       2 0.5;
       0 1;
       0.5 1;
       2 1]
icone=[1 2 5 4;
       2 3 6 5;
       4 5 8 7;
       5 6 9 8];
model.nnodes=9;
model.nelem=4;
model.kx=0;
model.ky=0;
model.c=10;
G=@(x,y) zeros(size(xnode,1)) + 50;
model.G=G(xnode(:,1),xnode(:,2));
model.ts=1;
model.rho=1;
model.cp=1;
model.maxit=5000;
model.tol=0.00000000001;
model.dt=0.01;
model.PHI_n=[0 0 0 0 0 0 0 0 0]';

DIR=[1 10;4 10;7 10];
NEU=[3 6 0;6 9 0];
ROB=[1 2 2 0;2 3 2 0;7 8 2 0;8 9 2 0];
PUN=[];

[PHI,Q]= fem2d_heat(xnode,icone,DIR,NEU,ROB,PUN,model);
fem2d_heat_graph_mesh(PHI,Q,xnode,icone,0,0);