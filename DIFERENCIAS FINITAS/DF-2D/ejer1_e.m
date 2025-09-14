#Ejercicio 1 - E
#Casuistica de diego en bordes donde chocan condiciones de borde:
#Robin se impone a neumann y dirichlet
#Dirichlet se impone a neumann
#dirichlet dirichlet = Promedio
#Neumann neuma -
clc;
clear all;
#Armado de la malla
x_size=3;
y_size=5;
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
model.rho=1;
model.cp=1;
model.maxit=100000;
model.tol=0.0000000000001;
model.k=2*ones(model.nnodes,1);
model.c=0*ones(model.nnodes,1);
model.ts=0; #FORWARD EULER
model.PHI_n=zeros(model.nnodes,1);
for j = 1:model.nnodes
      model.G(j)= 10*xnode(j,1)^2;
endfor

DIR=[1+ n*(n-1) 50;2+ n*(n-1) 50;3+ n*(n-1) 50;4+ n*(n-1) 50;5+ n*(n-1) 50;6+ n*(n-1) 50;7+ n*(n-1) 50; #Borde superior
     1+ n 0;1+ n*2 0;1+ n*3 0;1+ n*4 0;1+ n*5 0;         #Borde izquierdo
     1+ n*6 25; #Borde dir-dir = hacemos un promedio
     n*2  25;n*3  25;n*4  25;n*5  25;n*6  25; #Borde derecho
     n*n  75/2; #Borde dir-dir = hacemos un promedio
];
NEU=[];
ROB=[1 2 100 1;2 2 100 1;3 2 100 1;4 2 100 1;5 2 100 1;6 2 100 1;7 2 100 1];

[PHI,Q]=fdm2d(xnode,icone,DIR,NEU,ROB,model);
fdm2d_graph_mesh(PHI,Q,xnode,icone,0,0);
