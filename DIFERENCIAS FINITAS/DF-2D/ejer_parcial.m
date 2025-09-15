#Ejer_parcial
clc;
clear all;
model.nnodes=0;
#Armado de la malla
#Armado del rectangulo inferior                    % nodos por lado
xnode_x = -4:2:4;
xnode_y = -3:1:3;

[X,Y] = meshgrid(xnode_x,xnode_y);
X = X';
Y = Y';
xnode= [X(:), Y(:)] ;       % vector nodal

% --- generar icone ---
nelem = (length(xnode_x)-1)*(length(xnode_y)-1);          % cantidad de cuadrados
icone_inf = zeros(nelem,4);
k = 1;
for j = 1:(length(xnode_y)-1)               % filas (en y)
    for i = 1:(length(xnode_x)-1)           % columnas (en x)
        n1 = (j-1)*length(xnode_x) + i;     % nodo abajo-izq
        n2 = n1 + 1;          % nodo abajo-der
        n3 = n2 + length(xnode_x);          % nodo arriba-der
        n4 = n1 + length(xnode_x);          % nodo arriba-izq
        icone (k,:) = [n1 n2 n3 n4];
        k = k + 1;
    end
end

model.nnodes=size(xnode,1);
figure; hold on; axis equal;
patch('Faces',icone,'Vertices',xnode,'FaceColor','cyan','FaceAlpha',0.3);
plot(xnode(:,1), xnode(:,2), 'ko');  # nodos
indice_0=-1;
#Constantes del modelo :
model.k=3*ones(model.nnodes,1);
model.c=0*ones(model.nnodes,1);
model.PHI_n=zeros(model.nnodes,1);
for j = 1:model.nnodes
      model.G(j)= (100*xnode(j,1));
endfor
#Condiciones de borde
DIR=[];
ROB=[];
NEU=[];
for i=1:(model.nnodes)
  #Bordes faciles
  if(xnode(i,2)==-4)#Borde inferior
    NEU=[NEU;i 0 1];
  endif
  if(xnode(i,1)==3)#Borde Derecho
    ROB=[ROB;i 2 10 2];
  endif
  #Borde derecho
  if(xnode(i,2)==4)#Borde superior
    NEU=[NEU;i 0 3];
  endif
  #Borde superior
  if(xnode(i,1)==-4)
    DIR=[DIR;i 20];
  endif
  if(xnode(i,1)==0 && xnode(i,2)==0)
    indice_0=1;
  endif
endfor

model.ts=-1; #Estado estacionario
model.rho=0;
model.cp=0;
model.maxit=100000;
model.tol=0.0001;

[PHI,Q]=fdm2d(xnode,icone,DIR,NEU,ROB,model);
fdm2d_graph_mesh(PHI,Q,xnode,icone,0,0);




