function [localF] = fem2d_pstr_genF(nodes,Fg,th)
% Descripción: módulo para calcular el vector de fuerzas F para cada elemento,
% producto de la presencia de una fuerzas por unidad de volumen (generalmente 
% la gravedad) en dicho elemento. La integral se resuelve mediante cuadratura
% de punto medio, y se requiere evaluar el área del elemento.

% Entrada:
% * nodes: nodos (x,y) del elemento. Los elementos admisibles son de 3 o 4 nodos.
% * Fg: magnitud de la fuerza gravitatoria.
% * th: espesor de la placa.

% Salida:
% * localF: vector local de fuerzas.
% ----------------------------------------------------------------------
    if(size(nodes,1)==3)
        localF=zeros(6,1);
        localF(2)=1;
        localF(4)=1;
        localF(6)=1;
        A=0.5*abs(nodes(1,1)*(nodes(2,2)-nodes(3,2)) + nodes(2,1)*(nodes(3,2)-nodes(1,2)) + nodes(3,1)*(nodes(1,2)-nodes(2,2)));
        localF*=-A*th*Fg*(1/3)
    else
        localF=zeros(8,1);
        localF(2)=1;
        localF(4)=1;
        localF(6)=1;
        localF(8)=1;
        A=abs(nodes(1,1)*(nodes(2,2)-nodes(3,2)) + nodes(2,1)*(nodes(3,2)-nodes(1,2)) + nodes(3,1)*(nodes(1,2)-nodes(2,2)));
        localF*=-A*th*Fg*(1/4)
    end
end