#Ejercicio 1 - E
#Casuistica de diego en bordes donde chocan condiciones de borde:
#Robin se impone a neumann y dirichlet
#Dirichlet se impone a neumann
#dirichlet dirichlet = Promedio
#Neumann neumann -Promedio
clc;
clear all;
#Armado de la malla
x_size=2;
y_size=2;
n = 7;                        % nodos por lado
model.nnodes=n*n;
xv = linspace(0,x_size,n)
yv = linspace(0,y_size,n)

[X,Y] = meshgrid(xv,yv);
X = X';
Y = Y';
xnode = [X(:), Y(:)];         % vector nodal

% --- generar icone ---
nelem = (n-1)*(n-1);          % cantidad de cuadrados
icone = zeros(nelem,4);
k = 1;
for j = 1:(n-1)               % filas (en y)
    for i = 1:(n-1)           % columnas (en x)
        n1 = (j-1)*n + i;     % nodo abajo-izq
        n2 = n1 + 1;          % nodo abajo-der
        n3 = n2 + n;          % nodo arriba-der
        n4 = n1 + n;          % nodo arriba-izq
        icone(k,:) = [n1 n2 n3 n4];
        k = k + 1;
    end
end
#--------CONSTANTES----#
model.rho=2;
model.cp=1;
model.maxit=100000;
model.tol=1e-8;
model.k=2*ones(model.nnodes,1);
model.c=2*ones(model.nnodes,1);
model.ts=1; #FORWARD EULER
model.PHI_n=zeros(model.nnodes,1);
model.dt=0.25*min([abs(xnode(2,1)-xnode(1,1)) abs(xnode(2,2) - xnode(2,1))])/model.k(1);
for j = 1:model.nnodes
      model.G(j)= 0*xnode(j,1)^2+100;
endfor

DIR=[2 0;3 0;4 0;5 0;6 0;n 50; #Borde inferior
      n*2 100;n*3 100;n*4 100;n*5 100;n*6 100];
NEU=[];
ROB=[1 2 10 4;n+1 2 10 4;2*n+1 2 10 4;3*n+1 2 10 4;4*n+1 2 10 4;5*n+1 2 10 4; #Borde izquierdo
     6*n+1 2 10 4;#Esquina robin robin, donde coincide h y qinf pero no la normal. Se asume una normal
      6*n+2 2 10 3;6*n+3 2 10 3;6*n+4 2 10 3;6*n+4 2 10 3;6*n+5 2 10 3;6*n+6 2 10 3;n*n 2 10 3];#Borde Superior

[PHI,Q]=fdm2d(xnode,icone,DIR,NEU,ROB,model);
fdm2d_graph_mesh(PHI,Q,xnode,icone,0,0);
