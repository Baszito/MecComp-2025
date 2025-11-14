#Ejercicio 2- GUIA 2D
clc;
clear all;
xnode=[0 0;
       3/11 0;
       7/11 0;
       1 0;
       0 0.5;
       3/11 0.5;
       7/11 0.5;
       1 0.5;
       0 1;
       3/11 1;
       7/11 1;
       1 1];
       
icone = [
  1 2 6 5;   % E1
  2 3 7 6;   % E2
  3 4 8 7;   % E3
  5 6 10 9;  % E4
  6 7 11 10; % E5
  7 8 12 11; % E6
];


model.nnodes=12;
model.nelem=6;
model.kx=2;
model.ky=2;
model.c=0;
G=@(x,y) 10.*y;
model.G=G(xnode(:,1),xnode(:,2));
model.ts=-1;
model.rho=1;
model.cp=1;
model.maxit=1000;
model.tol=0.00000000001;
model.dt=0.005;
model.PHI_n=zeros(size(xnode,1),1);

% Dirichlet
DIR = [1 0;
       5 0;
       9 0];

% Robin
ROB = [4 8 10 5;
       8 12 10 5];

% Neumann
NEU = [9 10 0;
       10 11 0;
       11 12 0;
       1 2 1;
       2 3 1;
       3 4 1];

PUN=[];

[PHI,Q]= fem2d_heat(xnode,icone,DIR,NEU,ROB,PUN,model);
fem2d_heat_graph_mesh(PHI,Q,xnode,icone,0,4);