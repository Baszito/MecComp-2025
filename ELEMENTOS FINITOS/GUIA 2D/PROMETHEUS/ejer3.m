#Ejercicio 2- GUIA 2D
clc;
clear all;
xnode=[0 0;
       1 0;
       0.5 0.2667;
       0.5 0.8]
icone=[1 2 3 -1;
       1 3 4 -1;
       2 3 4 -1];

model.nnodes=4;
model.nelem=3;
model.kx=1;
model.ky=1;
model.c=10;
G=@(x,y) zeros(size(xnode,1));
model.G=G(xnode(:,1),xnode(:,2));
model.ts=-1;
model.rho=0;
model.cp=0;
model.maxit=1000;
model.tol=0.00001;
model.dt=1e-3;
model.PHI_n=[0 0 0 0 0 0 0 0 0]';

DIR=[];
NEU=[];
ROB=[1 2 5 0;1 4 5 0;2 4 5 0];
PUN=[];

[PHI,Q]= fem2d_heat(xnode,icone,DIR,NEU,ROB,PUN,model);
fem2d_heat_graph_mesh(PHI,Q,xnode,icone,0,0);