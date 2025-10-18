#Ejercicio SIMULACION - GUIA 2D
clc;
clear all;
xnode=[0 0;1/3 0;2/3 0;1 0;
       0 0.2;1/3 0.2;2/3 0.2;1 0.2];
icone=[1 2 6 5;2 3 7 6;3 4 8 7];
model.nnodes=8;
model.nelem=3;
model.kx=1;
model.ky=1;
model.c=0;
model.G=[100 100 100 100 100 100 100 100];
model.ts=-1;
model.rho=0;
model.cp=0;
model.maxit=1000;
model.tol=0.0001;
model.dt=0;
model.PHI_n=[0 0 0 0 0 0 0 0]';

DIR=[1 20;5 20;4 20;8 20];
NEU=[1 2 0;2 3 0;3 4 0];
PUN=[];
ROB=[5 6 5 0;6 7 5 0;7 8 5 0];

[PHI,Q]= fem2d_heat(xnode,icone,DIR,NEU,ROB,PUN,model);
fem2d_heat_graph_mesh(PHI,Q,xnode,icone,0,1);
axis equal
