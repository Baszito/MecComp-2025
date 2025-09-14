#Ejercicio 1 - D
#Casuistica de diego en bordes donde chocan condiciones de borde:
#Robin se impone a neumann y dirichlet
#Dirichlet se impone a neumann
#dirichlet dirichlet = Promedio
#Neumann neuma -
clc;
clear all;
#Armado de la malla
x_size=3;
y_size=3;
n = 6;                        % nodos por lado
model.nnodes=n*n;
xv = linspace(1,x_size,n)
yv = linspace(1,y_size,n)

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
model.rho=0;
model.cp=0;
model.maxit=1000;
model.tol=0.001;
model.k=1*ones(model.nnodes,1);
model.c=4*ones(model.nnodes,1);
model.ts=-1;
for j = 1:model.nnodes
      model.G(j)= 0*xnode(j,1);
endfor

DIR=[1 25;2 25;3 25;4 25;5 25;6 25];
NEU=[n*2 10 2;n*3 10 2;n*4 10 2;n*5 10 2;
     n+1 0 4;1+n*2 0 4;1+n*3 0 4;1+n*4 0 4;1+n*5 0 4
];
ROB=[1+n*5 2 50 3;2+n*5 2 50 3;3+n*5 2 50 3;4+n*5 2 50 3;5*n+5 2 50 3;n*n 2 50 3];

[PHI,Q]=fdm2d(xnode,icone,DIR,NEU,ROB,model);
fdm2d_graph_mesh(PHI,Q,xnode,icone,0,1);

