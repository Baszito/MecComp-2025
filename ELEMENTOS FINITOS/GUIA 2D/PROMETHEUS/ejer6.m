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
       
icone=[1 2 5 -1;
       2 3 5 -1;
       1 5 4 -1;
       3 6 5 -1;
       4 5 7 -1;
       5 6 9 -1;
       5 8 7 -1;
       5 9 8 -1];
model.nnodes=9;
model.nelem=8;
model.kx=1;
model.ky=1;
model.c=mean(100.*xnode(:,1));
G=@(x,y) 0.*xnode(:,1);
model.G=G(xnode(:,1),xnode(:,2));
model.ts=0;
model.rho=1;
model.cp=1;
model.maxit=1;
model.tol=0.00000000001;
model.dt=0.005;
model.PHI_n=sin(pi.*xnode(:,1)).*sin(pi.*xnode(:,2));

DIR=[1 0;2 0;3 0;4 0;6 0;7 0;8 0;9 0];
NEU=[];
ROB=[];
PUN=[];

[PHI,Q]= fem2d_heat(xnode,icone,DIR,NEU,ROB,PUN,model);
fem2d_heat_graph_mesh(PHI,Q,xnode,icone,0,1);