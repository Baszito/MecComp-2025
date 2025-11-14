clc;close all; clear all; more off;

xnode = [
0 0;
0 0.2;
0.333333333 0;
0.333333333 0.2;
0.66666666 0;
0.66666666 0.2;
1 0;
1 0.2;
];
xnode(:,2)=xnode(:,2)*2.5;

icone = [
1 3 2 -1;
3 4 2 -1;
3 5 6 4;
5 7 6 -1;
7 8 6 -1
];

Fixnodes = [
1 1 0;
1 2 0;
7 1 0;
7 2 0
];

Sideload = [
2 4 0 -1e4;
4 6 0 -1e4;
6 8 0 -1e4;
];

Pointload = [
2 0 -2000 0.2 0;
4 0 -2000 0.8 0;
];

disp('---------------------------------------------------------------');
disp('Inicializando modelo de datos...');

model.nnodes = size(xnode,1);
model.nelem = size(icone,1);

model.young = 200e6;
model.poiss = 0.3;
model.gravity = 0;
model.pstrs = 1;
model.thick = 0.02000000000000000;

disp('Iniciando el método numérico...');

% Llamada principal al Método de Elementos Finitos
[U,reaction,Def,Ten,Ten_VM] = fem2d_pstr(xnode,icone,Fixnodes,Sideload,Pointload,model);

disp('Finalizada la ejecución del método numérico.');

disp('---------------------------------------------------------------');
disp('Iniciando el post-procesamiento...');

% mode ---> modo de visualización:
%           [0] 2D - Con malla
%           [1] 2D - Sin malla
% graph --> tipo de gráfica:
%           [0] Estado inicial vs. desplazamiento
%           [1] Von Misses (escalar)
%           [2] Reacciones (vectorial)
% scale --> factor de escala para U (en veces)
mode = 0;
graph = 0;
scale = 1;
U = scale * U;
fem2d_pstr_graph_mesh(U,reaction,Ten_VM,xnode,icone,mode,graph);

disp('Finalizado el post-procesamiento.');
