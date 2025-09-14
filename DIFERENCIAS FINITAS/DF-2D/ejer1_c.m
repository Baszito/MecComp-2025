#Ejercicio 1 - C
clc;
clear all;
#Armado de la malla
x_size=4;
y_size=3;
n = 6;                        % nodos por lado
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
model.rho=0;
model.cp=0;
model.maxit=1000;
model.tol=0.001;
model.k=1*ones(model.nnodes,1);
model.c=0*ones(model.nnodes,1);
model.ts=-1;
for j = 1:model.nnodes
      model.G(j)= 2*xnode(j,1) ^2 + 3 *xnode(j,2)^3;
endfor

DIR=[1 50;n+1 50; 2*n+1 50; 3*n+1 50; 4*n+1 50; 5*n+1 50;
     5*n+2 50; 5*n+3 50;5*n+4 50;5*n+5 50;n*n 50
];
NEU=[1 0 1;2 0 1;3 0 1;4 0 1;5 0 1;6 0 1;
     n 0 2;n*2 0 2;n*3 0 2;n*4 0 2;n*5 0 2
];
ROB=[];

[PHI,Q]=fdm2d(xnode,icone,DIR,NEU,ROB,model);
fdm2d_graph_mesh(PHI,Q,xnode,icone,0,0);


