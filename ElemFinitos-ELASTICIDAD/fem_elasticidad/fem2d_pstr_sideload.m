function [F] = fem2d_pstr_sideload(F,Sideload,xnode,th)
% Descripción: módulo para calcular y ensamblar las contribuciones de pares de nodos
% (lados) donde se aplican cargas distribuidas.

% Entrada:
% * F: vector de fuerzas.
% * sideload: matriz con la información sobre fronteras con cargas distribuidas.
%   - Columnas 1-2: dos nodos contiguos formando un lado de un elemento.
%   - Columna 3: valor de fuerza en sentido eje-x.
%   - Columna 4: valor de fuerza en sentido eje-y.
% * xnode: matriz de nodos con pares (x,y) representando las coordenadas de cada nodo
%   de la malla.
% * th: espesor de la placa.

% Salida:
% * F: vector de fuerzas. Presenta modificaciones luego de aplicar la condición de borde.
% ----------------------------------------------------------------------
   for i=1:size(Sideload,1)
        %Longitud del lado
        L=sqrt(
        (xnode(Sideload(i,1),1) - xnode(Sideload(i,2),1))^2
        +(xnode(Sideload(i,1),2) - xnode(Sideload(i,2),2))^2);

        ind1=2*Sideload(i,1)-1;
        ind2=2*Sideload(i,2)-1;
        F(ind1)+=th*(Sideload(i,3)*L)/2;
        F(ind1+1)+=th*(Sideload(i,4)*L)/2;
        F(ind2)+=th*(Sideload(i,3)*L)/2;
        F(ind2+1)+=th*(Sideload(i,4)*L)/2;
    end
end
