function [localF] = fem2d_heat_genF(nodes,G)
% Descripción: módulo para calcular el vector de flujo térmico F para cada
% elemento, producto de la presencia de una fuente de calor en dicho elemento.
% La integral se resuelve mediante cuadratura de punto medio, y se requiere
% evaluar el área del elemento.

% Entrada:
% * nodes: nodos (x,y) del elemento. Los elementos admisibles son de 3 o 4 nodos.
% * G: fuente de calor.

% Salida:
% * localF: vector de flujo térmico (local).
% ----------------------------------------------------------------------
  n=size(nodes,1);
    if(n==3) %Elemento TRIANGULAR
        %Producto cruz
        v1=[nodes(2,1)-nodes(1,1) nodes(2,2)-nodes(1,2)];
        v2=[nodes(3,1)-nodes(1,1) nodes(3,2)-nodes(1,2)];
        area=abs(v1(1)*v2(2) - v1(2)*v2(1))/2;
        localF=[1;1;1];
        localF=(area*G/3)*localF;
    else
        %Producto cruz 1
        v1=[nodes(2,1)-nodes(1,1) nodes(2,2)-nodes(1,2)];
        v2=[nodes(3,1)-nodes(1,1) nodes(3,2)-nodes(1,2)];
        area1=abs(v1(1)*v2(2) - v1(2)*v2(1))/2;
        
        %Producto cruz 2
        v1=[nodes(4,1)-nodes(1,1) nodes(4,2)-nodes(1,2)];
        area2=abs(v1(1)*v2(2) - v1(2)*v2(1))/2;

        localF=[1;1;1;1];
        localF=((area1+area2)*G/4)*localF;
    end
end