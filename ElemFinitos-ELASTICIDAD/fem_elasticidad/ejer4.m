clc;close all; clear all; more off;

xnode = [
  0.0  0.0;
  0.2  0.0;
  0.4  0.0;
  0.0  0.2;
  0.2  0.2;
  0.4  0.2
];

icone = [
1 2 5 4;
2 3 6 5;
];

Fixnodes = [
1 1 0;1 2 0;4 1 0;4 2 0;2 2 0;
];

Sideload = [

];

Pointload = [
2   2000  0  0.1 0.4;
];

disp('---------------------------------------------------------------');
disp('Inicializando modelo de datos...');

model.nnodes = size(xnode,1);
model.nelem = size(icone,1);

model.young = 210e9;
model.poiss = 0.3;
model.gravity = -2.45e4;
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
