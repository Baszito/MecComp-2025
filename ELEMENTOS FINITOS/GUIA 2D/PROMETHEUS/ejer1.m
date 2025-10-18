#Ejercicio 1 - GUIA 2D

xnode=[0 0;1 0;0 1;0.5 0.5];
icone=[1 2 3 -1;1 3 4 -1];
model.nnodes=4;
model.nelem=2;
model.kx=1;
model.ky=1;
model.c=0;
model.G=[1 1];
model.ts=5;
model.rho=0;
model.cp=0;
model.maxit=1000;
model.tol=0.0001;
model.dt=0;
model.PHI_n=[0 0 0 0];

DIR=[4 0];
NEU=[1 3 5,1 2 10];
PUN=[];
ROB=[];

[PHI,Q]= fem2d_heat(xnode,icone,DIR,NEU,ROB,PUN,model);
fem2d_heat_graph_mesh(PHI,Q,xnode,icone,1,1);
