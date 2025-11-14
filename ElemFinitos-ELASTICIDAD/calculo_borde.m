% calculo Borde
% calculo_borde.m
% Lo probamos al carga_borde.m

clear; clc; close all;

% ------------------------------------------------------------------------
% Elegir caso: 'A' -> tracción constante, 'B' -> tracción variable
caseOption = 'B';
% ------------------------------------------------------------------------

% Definición del triángulo (nodos)
nodes = [0 0; 1 0; 0 1];   % [x1 y1; x2 y2; x3 y3]
edge = [1 2];              % lado entre nodo 1 y 2 (la base)

% Espesor (si corresponde) - si estás en tensión plana y querés incluir t:
t = 1.0; % por defecto 1; cambiá si querés efecto de espesor

switch upper(caseOption)
    case 'A'
        % A) tracción constante p = [0; -1000] N/m (por ejemplo)
        p_const = @(x,y) [0; -1000]; % N/m (lineal sobre el borde)
        fe = carga_borde(nodes, edge, p_const);
        fe = t * fe; % multiplicar por espesor si se usa (opcional)
        fprintf('CASO A: tracción constante p = [0; -1000] N/m\n');
        disp('Vector de fuerzas elementales fe (6x1):'); disp(fe);

        % comparación analítica: cada nodo del lado recibe L/2 * p
        L = sqrt( (nodes(edge(2),1)-nodes(edge(1),1))^2 + (nodes(edge(2),2)-nodes(edge(1),2))^2 );
        disp('Comparación analítica (cada nodo del lado = L/2 * p):');
        disp( (t * L/2) * [0; -1000] );

    case 'B'
        % B) tracción variable p = [0; -1000*x] por ejemplo
        p_var = @(x,y) [0; -1000 * x]; % N/m (depende de x)
        fe = carga_borde(nodes, edge, p_var);
        fe = t * fe; % multiplicar por espesor si se usa (segun sea TP o DP)
        fprintf('CASO B: tracción variable p = [0; -1000*x] N/m\n');
        disp('Vector de fuerzas elementales fe (6x1):'); disp(fe);

    otherwise
        error('caseOption debe ser ''A'' o ''B''');
end
