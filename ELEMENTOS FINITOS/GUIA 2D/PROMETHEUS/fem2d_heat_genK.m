function [localK] = fem2d_heat_genK(nodes,kx,ky)
% Descripción: módulo para calcular y evaluar de forma numérica la matriz de
% difusión K. Se utilizan funciones de forma en coordenadas naturales y se 
% resuelve la integral de forma numérica utilizando cuadratura de Gauss.

% Entrada:
% * nodes: nodos (x,y) del elemento. Los elementos admisibles son de 3 o 4 nodos.
% * kx: conductividad térmica orientada en eje-x.
% * ky: conductividad térmica orientada en eje-y.

% Salida:
% * localK: matriz de difusión del elemento (local).
% ----------------------------------------------------------------------
    n=size(nodes,1); 
    Kdif=[kx 0;0 ky];
    if(n==3)%aca ta todo joya no hay que integrar
        J=[nodes(2,1)-nodes(1,1) nodes(2,2)-nodes(1,2);
            nodes(3,1)-nodes(1,1) nodes(3,2)-nodes(1,2)];
        DN=[1 -1 0;1 0 -1];
        
        B=inv(J)*DN;
        localK=det(J)*(0.5)*B'*Kdif*B;
    else
        %usamos cuadratura de gauss-legendre con dos puntos en x e y, 4 totales
        factor=1/sqrt(3);
        inter=[-factor factor];
        localK=zeros(4,4);
        for i=1:2
            for j=1:2
                s=inter(i);
                t=inter(j);
                DN=[-(1-t) (1-t) (1+t) -(1+t);
                    -(1-s) -(1+s) (1+s) (1-s)];
                DN=1/4 * DN;
                J=DN*nodes;
                B=inv(J)*DN;
                localK+= B'*Kdif*B*det(J);
            end
        end
    end
end