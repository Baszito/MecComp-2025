close all; clear all; more off;

xnode = [
  0.0000000000000000, 0.0000000000000000;
  0, 0.1000000000000000;
  1, 0.0000000000000000;
  1, 0.1000000000000000
];

icone = [
       1,      3,      4,     2;
];

Fixnodes = [
       1,      1, 0.0000000000000000;
       1,      2, 0.0000000000000000;
       2,      1, 0.0000000000000000;
       2,      2, 0.0000000000000000;
];

Sideload = [
];

Pointload = [
      1,0,-1000,1,0.05
];

disp('---------------------------------------------------------------');
disp('Inicializando modelo de datos...');

model.nnodes = size(xnode,1);
model.nelem = size(icone,1);

model.young = 21000000000000.0000000000000000;
model.poiss = 0.3000000000000000;
model.gravity = 9.8100000000000005;
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
mode = 1;
graph = 0;
scale = 1;
U = scale * U;
fem2d_pstr_graph_mesh(U,reaction,Ten_VM,xnode,icone,mode,graph);

disp('Finalizado el post-procesamiento.');
