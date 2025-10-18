clc;
clear all;

%% =================== DATOS ===================
he = 1/3;        % longitud de cada elemento
k = 2;           % conductividad / módulo
n_nodes = 4;     % nodos globales
n_elem = 3;      % número de elementos

%% =================== FUNCIONES DE FORMA ===================
% Funciones de forma cuadráticas en coordenadas naturales xi [-1,1]
N1 = @(xi) xi.*(xi-1)/2;
N2 = @(xi) 1 - xi.^2;
N3 = @(xi) xi.*(xi+1)/2;

% Derivadas respecto a xi
dN1 = @(xi) xi - 1/2;
dN2 = @(xi) -2*xi;
dN3 = @(xi) xi + 1/2;

% Derivadas respecto a x
dN1dx = @(xi) (2/he) * dN1(xi);
dN2dx = @(xi) (2/he) * dN2(xi);
dN3dx = @(xi) (2/he) * dN3(xi);

dN = {dN1dx, dN2dx, dN3dx};

%% =================== MAPEO DE NODOS ===================
% Cada fila: nodos globales del elemento
element_nodes = [1 2 3;
                 2 3 4;
                 3 4 4];  % último nodo repetido para simplificar

%% =================== MATRIZ GLOBAL ===================
K = zeros(n_nodes,n_nodes);
F = zeros(n_nodes,1);   % vector fuente, si hay (puede ser cero)

%% =================== ENSAMBLAJE ===================
for e = 1:n_elem
    nodes = element_nodes(e,:);

    % MATRIZ LOCAL 3x3
    Ke = zeros(3,3);
    for i = 1:3
        for j = 1:3
            % Integral con cuadratura de Gauss 2 puntos
            Ke(i,j) = k * cuad_gauss_c(@(xi) dN{i}(xi).*dN{j}(xi), -1, 1, 1, 2);
        end
    end

    % ENSAMBLAJE EN K GLOBAL
    for i = 1:3
        for j = 1:3
            g_i = nodes(i);
            g_j = nodes(j);
            if g_i <= n_nodes && g_j <= n_nodes
                K(g_i,g_j) = K(g_i,g_j) + Ke(i,j);
            end
        end
    end

    % VECTOR FUENTE LOCAL (si hay alguna fuente f(x), ejemplo f=1)
    fe_local = zeros(3,1);
    for i = 1:3
        fe_local(i) = cuad_gauss_c(@(xi) N1(xi), -1, 1, 1, 2);  % ejemplo: usando N1, cambiar según tu fuente
    end

    % ENSAMBLAJE EN F GLOBAL
    for i = 1:3
        g_i = nodes(i);
        if g_i <= n_nodes
            F(g_i) = F(g_i) + fe_local(i);
        end
    end
end

%% =================== RESULTADOS ===================
disp('Matriz de rigidez K:');
disp(K);
disp('Vector de fuerzas F:');
disp(F);

