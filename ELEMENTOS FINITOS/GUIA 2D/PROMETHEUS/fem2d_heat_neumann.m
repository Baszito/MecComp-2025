function [F] = fem2d_heat_neumann(F,NEU,xnode)
% Descripción: módulo para calcular y ensamblar las contribuciones de nodos
% pertenecientes a fronteras de tipo Neumann.

% Entrada:
% * F: vector de flujo térmico.
% * NEU: matriz con la información sobre la frontera de tipo Neumann. 
%   - Columnas 1-2: dos nodos contiguos formando un lado de un elemento.
%   - Columna 3: valor de flujo térmico (q) asociado al lado del elemento.
% * xnode: matriz de nodos con pares (x,y) representando las coordenadas de 
%   cada nodo de la malla.

% Salida:
% * F: vector de flujo térmico. Presenta modificaciones luego de aplicar la condición de borde.
% --------------------------------------------------------------------
  n = size(NEU, 1);

  for i = 1:n
    n1 = NEU(i,1);
    n2 = NEU(i,2);
    q  = NEU(i,3);

    % Longitud del lado
    L = sqrt((xnode(n2,1)-xnode(n1,1))^2 + (xnode(n2,2)-xnode(n1,2))^2);

    % Ensamblaje (aportan ambos nodos)
    F(n1) -= (q * L) / 2;
    F(n2) -= (q * L) / 2;
  end
end