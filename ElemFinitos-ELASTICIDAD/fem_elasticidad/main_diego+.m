clc;close all; clear all; more off;

xnode = [
  0.0000000000000000, 0.0000000000000000;
  0.0000000000000000, 0.2000000000000000;
  0.3333333333333333, 0.0000000000000000;
  0.3333333333333333, 0.2000000000000000;
  0.6666666666666666, 0.0000000000000000;
  0.6666666666666666, 0.2000000000000000;
  1.0000000000000000, 0.0000000000000000;
  1.0000000000000000, 0.2000000000000000;
];

xnode(:,2) = xnode(:,2)*2.5;

icone = [
       1,      3,      2,     -1;
       3,      4,      2,     -1;
       3,      5,      6,      4;
       5,      7,      6,     -1;
       7,      8,      6,     -1;
];

Fixnodes = [
       1,      1, 0.0000000000000000;
       1,      2, 0.0000000000000000;
       7,      1, 0.0000000000000000;
       7,      2, 0.0000000000000000;
];

Sideload = [
       4,      2, 0.0000000000000000, -10000.0000000000000000;
       6,      4, 0.0000000000000000, -10000.0000000000000000;
       8,      6, 0.0000000000000000, -10000.0000000000000000;
];

Pointload = [
       1, 0.0000000000000000, -2000.0000000000000000, 0.2000000000000000, 0.0000000000000000
       4, 0.0000000000000000, -2000.0000000000000000, 0.8000000000000000, 0.0000000000000000
];

disp('---------------------------------------------------------------');
disp('Inicializando modelo de datos...');

model.nnodes = size(xnode,1);
model.nelem = size(icone,1);

model.young = 200000000;
model.poiss = 0.3000000000000000;
model.gravity = 0.0000000000000000;
model.pstrs = 1;
model.thick = 0.0200000000000000;

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
