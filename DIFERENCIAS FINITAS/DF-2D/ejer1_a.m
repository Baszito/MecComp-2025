##GUIA 2 - Armado de modelos para comprobar
clc;
clear all;
#-------------------------------1-------------------------------#

#EJERCICIO A)
model.nnodes=4;
model.rho=0;
model.cp=0;
model.maxit=1000;
model.tol=0.001;
model.k=2*ones(model.nnodes,1);
model.c=0*ones(model.nnodes,1);
model.G=300*ones(model.nnodes,1);

#Problema estacionario
model.ts=-1;#Estacionario
DIR=[1 20;2 20;3 20; 4 20];
NEU=[];
ROB=[];

#Primera malla, solo 4 nodos, uniforme
xnode=[0 0;1 0;1 1; 0 1];
icone=[1 2 3 4];
model.nnodes=4;

[PHI,Q]=fdm2d(xnode,icone,DIR,NEU,ROB,model);
fdm2d_graph_mesh(PHI,Q,xnode,icone,0,0);

#Segunda malla, 9 nodos uniforme
model.nnodes=9;
model.rho=0;
model.cp=0;
model.maxit=1000;
model.tol=0.001;
model.k=2*ones(model.nnodes,1);
model.c=0*ones(model.nnodes,1);
model.G=300*ones(model.nnodes,1);

#Problema estacionario
model.ts=-1;#Estacionario
DIR=[1 20;2 20;3 20; 4 20;6 20;7 20;8 20;9 20];
NEU=[];
ROB=[];

xnode=[0 0;0.5 0;1 0;0 0.5;0.5 0.5;1 0.5;0 1; 0.5 1;1 1];
icone=[1 2 5 4;2 3 6 5;4 5 8 7;5 6 9 8];

[PHI,Q]=fdm2d(xnode,icone,DIR,NEU,ROB,model);
fdm2d_graph_mesh(PHI,Q,xnode,icone,0,0);

#Tercera malla, 9 nodos no-uniforme
model.nnodes=9;
model.rho=0;
model.cp=0;
model.maxit=1000;
model.tol=0.001;
model.k=2*ones(model.nnodes,1);
model.c=0*ones(model.nnodes,1);
model.G=300*ones(model.nnodes,1);

#Problema estacionario
model.ts=-1;#Estacionario
DIR=[1 20;2 20;3 20; 4 20;6 20;7 20;8 20;9 20];
NEU=[];
ROB=[];

xnode=[0 0;1/3 0;1 0;0 2/3;1/3 2/3;1 2/3;0 1; 1/3 1;1 1];
icone=[1 2 5 4;2 3 6 5;4 5 8 7;5 6 9 8];

[PHI,Q]=fdm2d(xnode,icone,DIR,NEU,ROB,model);
fdm2d_graph_mesh(PHI,Q,xnode,icone,0,0);

#Cuarta malla, 36 nodos uniformes
model.nnodes=36;
model.rho=0;
model.cp=0;
model.maxit=1000;
model.tol=0.001;
model.k=2*ones(model.nnodes,1);
model.c=0*ones(model.nnodes,1);
model.G=300*ones(model.nnodes,1);
s=6;
#Problema estacionario
model.ts=-1;#Estacionario
[X,Y]=meshgrid(linspace(0,1,s),linspace(0,1,s));
X=X';Y=Y';
xnode=[X(:),Y(:)];
DIR=[1 20; 2 20; 3 20; 4 20; 5 20; s 20;
     2*s 20; 3*s 20; 4*s 20; 5*s 20; 6*s 20;
     s+1 20; 2*s+1 20; 3*s+1 20; 4*s+1 20; 5*s+1 20;
     5*s+2 20; 5*s+3 20; 5*s+4 20; 5*s+5 20]
NEU=[];
ROB=[];


icone = [
 1  2  8  7;
 2  3  9  8;
 3  4 10  9;
 4  5 11 10;
 5  6 12 11;
 7  8 14 13;
 8  9 15 14;
 9 10 16 15;
10 11 17 16;
11 12 18 17;
13 14 20 19;
14 15 21 20;
15 16 22 21;
16 17 23 22;
17 18 24 23;
19 20 26 25;
20 21 27 26;
21 22 28 27;
22 23 29 28;
23 24 30 29;
25 26 32 31;
26 27 33 32;
27 28 34 33;
28 29 35 34;
29 30 36 35
];

[PHI,Q]=fdm2d(xnode,icone,DIR,NEU,ROB,model);
fdm2d_graph_mesh(PHI,Q,xnode,icone,0,0);


#EJERCICIO B)
#model.nnodes=4;
#model.rho=0;
#model.cp=0;
#model.maxit=1000;
#model.tol=0.001;
#model.k=2*ones(model.nnodes,1);
#model.c=4*ones(model.nnodes,1);
#model.G=300*ones(model.nnodes,1);

#Primera malla, solo 4 nodos, uniforme
#xnode=[0 0;1 0;1 2; 0 2];
#icone=[1 2 3 4];
#model.nnodes=4;

#Problema estacionario
#model.ts=-1;#Estacionario
#DIR=[1 20;2 0;3 20; 4 0];
#NEU=[];
#ROB=[];
#[PHI,Q]=fdm2d(xnode,icone,DIR,NEU,ROB,model);
#fdm2d_graph_mesh(PHI,Q,xnode,icone,0,0);



