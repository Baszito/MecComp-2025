function [K,F] = fem2d_heat_dirichlet(K,F,DIR)
% Descripción: módulo para calcular y ensamblar las contribuciones de nodos
% pertenecientes a fronteras de tipo Dirichlet.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * DIR: matriz con la información sobre la frontera de tipo Dirchlet.
%   - Columna 1: número de nodo.
%   - Columna 2: valor en ese nodo (escalar).

% Salida:
% * K: matriz del sistema (difusión + reacción) con modificaciones luego de
%   aplicar la condición de borde.
% * F: vector de flujo térmico con modificaciones luego de aplicar la condición
%   de borde.
% ----------------------------------------------------------------------
    n=length(DIR(:,1));
    for i=1:n
       K(DIR(i,1),:)=zeros(1,length(K(1,:)));
       K(DIR(i,1),DIR(i,1))=1;
       F(DIR(i,1))=DIR(i,2);
    end
end