function [localK] = fem2d_pstr_genK(nodes,D,th)
% Descripción: módulo para calcular y evaluar de forma numérica la matriz de rigidez K.
% Se utilizan funciones de forma en coordenadas naturales y se resuelve la integral
% de forma numérica utilizando cuadratura de Gauss.

% Entrada:
% * nodes: nodos (x,y) del elemento. Los elementos admisibles son de 3 o 4 nodos.
% * D: matriz constitutiva.
% * t: espesor de la placa.

% Salida:
% * localK: matriz de rigidez del elemento (local).
% ----------------------------------------------------------------------
  %Acordate en caso de tension plana multiplciar por el espesor
%---------------CASO TRIANGULITO---------------%
  %Area del triangulo, lo usamos en todas las B
  if(size(nodes,1)==3)
      A=0.5*abs(nodes(1,1)*(nodes(2,2)-nodes(3,2)) + nodes(2,1)*(nodes(3,2)-nodes(1,2)) + nodes(3,1)*(nodes(1,2)-nodes(2,2)));
      b(1)=nodes(2,2)-nodes(3,2);
      b(2)=nodes(3,2)-nodes(1,2);
      b(3)=nodes(1,2)-nodes(2,2);

      c(1)=nodes(3,1)-nodes(2,1);
      c(2)=nodes(1,1)-nodes(3,1);
      c(3)=nodes(2,1)-nodes(1,1);

      B1=zeros(3,2);
      B1=[b(1) 0;0 c(1);c(1) b(1)]
      B2=zeros(3,2);
      B2=[b(2) 0;0 c(2);c(2) b(2)]
      B3=zeros(3,2);
      B3=[b(3) 0;0 c(3);c(3) b(3)]

      B=[B1 B2 B3];

      localK=(1/(4*A)).*(B' * D * B);
   else
      %Elemento cuadrangular
      %Integracion, con 4 puntitos de gauss,w=1
        factor=1/sqrt(3);
        inter=[-factor factor];
        localK=zeros(8,8);
        for i=1:2
            for j=1:2
                s=inter(i);
                t=inter(j);

                DN = [ -(1-t)  (1-t) (1+t) -(1+t);
                 -(1-s) -(1+s) (1+s)  (1-s) ] / 4;

              J = DN * nodes;

              dN = inv(J) * DN;

              B = zeros(3,8);
              B(1,1:2:8) = dN(1,:);
              B(2,2:2:8) = dN(2,:);
              B(3,1:2:8) = dN(2,:);
              B(3,2:2:8) = dN(1,:);

              localK += B' * D * B * det(J);
            end
        end

   end
    localK*=th;
end
