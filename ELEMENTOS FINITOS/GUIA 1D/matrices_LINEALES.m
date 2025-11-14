clc;
clear all;

%% =================== DATOS ===================
h = 2/3;         % longitud de cada elemento
k = 1;           % conductividad / módulo
c = 1;           % coeficiente de reacción
G = 0;           % fuente (si la hay)
n_nodes = 4;     % cantidad de nodos globales
n_elem = 3;      % cantidad de elementos

%% =================== MAPEO DE NODOS ===================
% Cada fila indica los nodos globales que forman el elemento e
element_nodes = [1 2;
                 2 3;
                 3 4];

%% =================== MATRICES GLOBALES ===================
K = zeros(n_nodes,n_nodes);
C = zeros(n_nodes,n_nodes);
F = zeros(n_nodes,1);

%% =================== ENSAMBLAJE ===================
for e = 1:n_elem
    nodes = element_nodes(e,:);

    % --- MATRICES LOCALES ---
    Ke = (k/h) * [ 1 -1; -1 1 ];
    Ce = (c*h/6) * [ 2 1; 1 2 ];
    Fe = (G*h/2) * [ 1; 1 ];

    % --- ENSAMBLAJE GLOBAL ---
    for i = 1:2
        for j = 1:2
            g_i = nodes(i);
            g_j = nodes(j);
            K(g_i,g_j) = K(g_i,g_j) + Ke(i,j);
            C(g_i,g_j) = C(g_i,g_j) + Ce(i,j);
        end
    end

    % --- VECTOR FUENTE ---
    for i = 1:2
        g_i = nodes(i);
        F(g_i) = F(g_i) + Fe(i);
    end
end

%% =================== RESULTADOS ===================
disp('Matriz de rigidez K:');
disp(K);
disp('Matriz de reacción C:');
disp(C);
disp('Vector de fuerzas F:');
disp(F);

%% (Opcional) Sistema global total:
% A = K + C;
% T = A\F;
% disp('Temperaturas nodales:');
% disp(T);

