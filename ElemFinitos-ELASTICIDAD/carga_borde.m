function fe = carga_borde(nodes, edge, p_var)
% Calcula la carga nodal equivalente (6x1) de un triángulo lineal (CST)
% debida a una tracción vectorial p(x,y) aplicada sobre el lado 'edge'.
%
% nodes : 3x2  => [x1 y1; x2 y2; x3 y3]
% edge  : 1x2  => índices de nodos que definen el lado, por ej. [1 2]
% p_var : function var @(x,y) -> [px; py]
%
% fe : 6x1  vector de fuerzas elementales (orden nodal [u1 v1 u2 v2 u3 v3]')

% puntos de Gauss para [0,1] (2-puntos)
gp = [0.211324865405187, 0.788675134594813];
gw = [0.5, 0.5];

% índice del tercer nodo
allIdx = 1:3;
k = setdiff(allIdx, edge);

% coordenadas de nodos del lado
xi = nodes(edge(1),1); yi = nodes(edge(1),2);
xj = nodes(edge(2),1); yj = nodes(edge(2),2);

% longitud del lado
L = sqrt( (xj-xi)^2 + (yj-yi)^2 );

% inicializar fe local del borde (4x1: [uxi; uyi; uxj; uyj])
fe_edge = zeros(4,1);

for q = 1:length(gp)
    s = gp(q);           % parámetro en [0,1]
    % punto físico en el lado
    x_s = xi + s*(xj - xi);
    y_s = yi + s*(yj - yi);

    % funciones de forma en el lado (lineales)
    N_i = 1 - s;
    N_j = s;

    % tracción en (x_s, y_s)
    p = p_var(x_s, y_s); % debe devolver [px; py] (2x1)
    % Control
    if numel(p)~=2
        error('p_var debe devolver 2 componentes [px;py]');
    end

    % CUADRATURA
    % contribución al integrando -> N(s) * p(x_s) * |J|, J = L (ds/dparam)
    w = gw(q);
    % sumar: para nodo i
    fe_edge(1:2) = fe_edge(1:2) + N_i * p * w * L;
    % para nodo j
    fe_edge(3:4) = fe_edge(3:4) + N_j * p * w * L;
end

% ahora armar el vector fe (6x1) del elemento triángulo CST
fe = zeros(6,1);
% nodo i -> posiciones 1,2
fe( (edge(1)-1)*2 + (1:2) ) = fe( (edge(1)-1)*2 + (1:2) ) + fe_edge(1:2);
% nodo j -> posiciones 1,2 offset por j
fe( (edge(2)-1)*2 + (1:2) ) = fe( (edge(2)-1)*2 + (1:2) ) + fe_edge(3:4);

% nodo k (tercer nodo) queda con ceros si no hay tracción en sus bordes
end

