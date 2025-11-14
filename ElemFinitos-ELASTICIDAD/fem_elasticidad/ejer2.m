clc;close all; clear all; more off;

xnode = [
  0.0  0.0;
  0.6  0.0;
  0.25 0.15;
  0.35 0.15;
  0.25 0.25;
  0.35 0.25;
  0.0  0.4;
  0.6  0.4;
];

icone = [
1 2 4 3;
5 6 8 7;
4 2 8 6;
1 3 5 7;
];

Fixnodes = [
1 1 0;1 2 0;7 1 0;7 2 0
];

Sideload = [
  2   8  0.5e6  0;
];

Pointload = [];

disp('---------------------------------------------------------------');
disp('Inicializando modelo de datos...');

model.nnodes = size(xnode,1);
model.nelem = size(icone,1);

model.young = 70e9;
model.poiss = 0.33;
model.gravity = 0;
model.pstrs = 1;
model.thick = 0.01000000000000000;

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

U
Ten
disp("Concentracion de tensiones en las aristas")
Ten_VM./0.5e6
