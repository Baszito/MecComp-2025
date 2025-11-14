#Ejercicio 8- GUIA 2D
clc;
clear all;
xnode=[-1 -1;
       1 -1;
       -1.5 0;
       1.5 0;
       -1 1;
       1 1];
       
icone = [
  1 2 6 5;   % E1
  1 5 3 -1;   % E2
  2 4 6 -1;   % E3
];


model.nnodes=6;
model.nelem=3;
model.kx=1;
model.ky=1;
model.c=0;
G=@(x,y) 50+0.*(x);
model.G=G(xnode(:,1),xnode(:,2));
model.ts=1;
model.rho=1;
model.cp=1;
model.maxit=100;
model.tol=0.00000000001;
model.dt=0.01;
model.PHI_n=zeros(size(xnode,1),1);

DIR=[1 0;2 0];
NEU=[1 3 0;2 4 0;3 5 0;4 6 0];
ROB=[];
PUN=[];

[PHI,Q]= fem2d_heat(xnode,icone,DIR,NEU,ROB,PUN,model);
fem2d_heat_graph_mesh(PHI,Q,xnode,icone,0,0);