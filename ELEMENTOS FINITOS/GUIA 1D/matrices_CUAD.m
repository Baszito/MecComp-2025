clc;
clear all;

%% =================== DATOS ===================
he = 2/3;        % longitud de cada elemento
k = 1;           % conductividad / módulo
c = 1;           % coeficiente de reacción
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
N  = {N1, N2, N3};

%% =================== MAPEO DE NODOS ===================
element_nodes = [1 2 3;
                 2 3 4;
                 3 4 4];  % último nodo repetido para simplificar

%% =================== MATRICES GLOBALES ===================
K = zeros(n_nodes,n_nodes);
C = zeros(n_nodes,n_nodes);
F = zeros(n_nodes,1);

%% =================== ENSAMBLAJE ===================
for e = 1:n_elem
    nodes = element_nodes(e,:);

    % MATRICES LOCALES
    Ke = zeros(3,3);
    Ce = zeros(3,3);

    % --- MATRIZ DE RIGIDEZ (conducción) ---
    for i = 1:3
        for j = 1:3
            Ke(i,j) = k * cuad_gauss_c(@(xi) dN{i}(xi).*dN{j}(xi), -1, 1, 1, 2);
        end
    end

    % --- MATRIZ DE REACCIÓN (masa) ---
    for i = 1:3
        for j = 1:3
            Ce(i,j) = c * (he/2) * cuad_gauss_c(@(xi) N{i}(xi).*N{j}(xi), -1, 1, 1, 2);
        end
    end

    % --- ENSAMBLAJE GLOBAL ---
    for i = 1:3
        for j = 1:3
            g_i = nodes(i);
            g_j = nodes(j);
            if g_i <= n_nodes && g_j <= n_nodes
                K(g_i,g_j) = K(g_i,g_j) + Ke(i,j);
                C(g_i,g_j) = C(g_i,g_j) + Ce(i,j);
            end
        end
    end

    % --- VECTOR FUENTE LOCAL (ejemplo f=1) ---
    fe_local = zeros(3,1);
    for i = 1:3
        fe_local(i) = (he/2) * cuad_gauss_c(@(xi) N{i}(xi), -1, 1, 1, 2);
    end

    % --- ENSAMBLAJE GLOBAL DE F ---
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
disp('Matriz de reacción C:');
disp(C);
disp('Vector de fuerzas F:');
disp(F);


