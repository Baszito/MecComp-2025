#Ejercicio 2 - D
clc;
clear all;
model.nnodes=0;
#Armado de la malla
#Armado del rectangulo inferior                    % nodos por lado
xnode_inf_x = 0 : 3/48 : 3;
xnode_inf_y = 0 : 2/14 : 2;

[X,Y] = meshgrid(xnode_inf_x,xnode_inf_y);
X = X';
Y = Y';
xnode_inf = [X(:), Y(:)] ;       % vector nodal

% --- generar icone ---
nelem = (length(xnode_inf_x)-1)*(length(xnode_inf_y )-1);          % cantidad de cuadrados
icone_inf = zeros(nelem,4);
k = 1;
for j = 1:(length(xnode_inf_y )-1)               % filas (en y)
    for i = 1:(length(xnode_inf_x)-1)           % columnas (en x)
        n1 = (j-1)*length(xnode_inf_x) + i;     % nodo abajo-izq
        n2 = n1 + 1;          % nodo abajo-der
        n3 = n2 + length(xnode_inf_x);          % nodo arriba-der
        n4 = n1 + length(xnode_inf_x);          % nodo arriba-izq
        icone_inf (k,:) = [n1 n2 n3 n4];
        k = k + 1;
    end
end

#Armado del rectangulo superior                    % nodos por lado
xnode_sup_x = 0 : 1/16 : 1;
xnode_sup_y = 2 + (1:16)/8;
[X,Y] = meshgrid(xnode_sup_x,xnode_sup_y);
X = X';
Y = Y';
xnode_sup = [X(:), Y(:)] ;       % vector nodal

#Armado de la malla final
xnode=[xnode_inf;xnode_sup];
model.nnodes = size(xnode,1);

% --- generar icone ---
nelem = (length(xnode_sup_x)-1)*(length(xnode_sup_y)-1);          % cantidad de cuadrados
icone_sup = zeros(nelem,4);
k = 1;
for j = 1:(length(xnode_sup_y)-1)               % filas (en y)
    for i = 1:(length(xnode_sup_x)-1)           % columnas (en x)
        n1 = (j-1)*length(xnode_sup_x) + i;     % nodo abajo-izq
        n2 = n1 + 1;          % nodo abajo-der
        n3 = n2 + length(xnode_sup_x);          % nodo arriba-der
        n4 = n1 + length(xnode_sup_x);          % nodo arriba-izq
        icone_sup(k,:) = [n1 n2 n3 n4];
        k = k + 1;
    end
end

#Armado final de icone
icone_sup=icone_sup + size(xnode_inf,1); #corremos los indices
#Solo faltan los nodos que van entre ambas figuras :
nelem = (length(xnode_sup_x));         % cantidad de cuadrados
icone_mid = zeros(nelem-1,4);
k = 1;
for i = 1:(nelem-1)           % columnas (en x)
        n1 = (size(xnode_inf,1)-length(xnode_inf_x))+i;     % nodo abajo-izq
        n2 = n1 + 1;          % nodo abajo-der
        n3 = n2 + (length(xnode_inf_x));          % nodo arriba-der
        n4 = n1 + (length(xnode_inf_x));          % nodo arriba-izq
        icone_mid(k,:) = [n1 n2 n3 n4];
        k = k + 1;
end
icone=[icone_inf ;icone_mid; icone_sup];
#Grafiquita, para ver que andemos bien

figure; hold on; axis equal;
patch('Faces',icone,'Vertices',xnode,'FaceColor','cyan','FaceAlpha',0.3);
plot(xnode(:,1), xnode(:,2), 'ko');  # nodos

#Constantes del modelo :
model.k=1*ones(model.nnodes,1);
model.c=0*ones(model.nnodes,1);
model.PHI_n=zeros(model.nnodes,1);
for j = 1:model.nnodes
      model.G(j)= (0*xnode(j,1))+100;
endfor

#Condiciones de borde
DIR=[];
ROB=[];
NEU=[];
for i=1:(model.nnodes)
  #Bordes faciles
  if(xnode(i,2)==0)#Borde inferior
    ROB=[ROB;i 2 50 1];
  endif
  if(xnode(i,1)==0)#Borde IZQUIERDO
    DIR=[DIR;i 20];
  endif
  #Borde derecho corto
  if(xnode(i,1)==3)
    NEU=[NEU;i 0 2];
  endif
  #Borde derecho largo
  if(xnode(i,1)==1 && xnode(i,2)>2)
    NEU=[NEU;i 0 2];
  endif
  #Borde superior
  if(xnode(i,2)==6)
    DIR=[DIR;i 0];
  endif
endfor
#ITEM A) - ESTADO ESTACIONARIO
model.ts=-1; #Estado estacionario
model.rho=0;
model.cp=0;
model.maxit=100000;
model.tol=0.0001;

[PHI,Q]=fdm2d(xnode,icone,DIR,NEU,ROB,model);
fdm2d_graph_mesh(PHI,Q,xnode,icone,0,0);

#ITEM D) - ESTADO TRANSITORIO
#for j = 1:model.nnodes
#      model.G(j)= (0*xnode(j,1))+100;
#endfor

#model.ts=1; #BACKWARD EULER HASTA 2 SEG
#model.rho=1;
#model.cp=1;
#model.dt=0.05;
#model.maxit=1000;
#model.tol=1e-5;

#[PHI,Q]=fdm2d(xnode,icone,DIR,NEU,ROB,model); #Uso este PHI Y Q
#fdm2d_graph_mesh(PHI,Q,xnode,icone,0,0);
