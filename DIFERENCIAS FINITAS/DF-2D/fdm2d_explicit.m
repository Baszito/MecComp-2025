function [PHI, Q] = fdm2d_explicit(K,F,xnode,neighb,model,dt)
% Descripción: módulo para resolver el sistema lineal de ecuaciones utilizando esquema
% temporal implícito.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * xnode: matriz de nodos con pares (x,y) representando las coordenadas de cada nodo
% de la malla.
% * neighb: matriz de vecindad.
% * model: struct con todos los datos del modelo (constantes, esquema numérico, etc.)
% * dt: paso temporal crítico para método explícito.

% Salida:
% * PHI: matriz solución. Cada elemento del vector representa un valor escalar
%   asociado a cada nodo de la malla, y su posición dentro del vector depende de
%   cómo se especificó cada nodo en xnode. Cada columna representa una iteración
%   del esquema temporal (en total nit columnas).
% * Q: matriz de flujo de calor. Para cada nodo se halla un vector bidimensional
%   de flujo de calor, representado por un par (Qx,Qy). Cada par de columnas
%   representa una iteración del esquema temporal (en total 2×nit columnas).
% ----------------------------------------------------------------------

    A = dt/(model.rho*model.cp);

    I = eye(model.nnodes,model.nnodes);

    PHI = model.PHI_n;
    PHI_n = model.PHI_n;
    Q = fdm2d_flux(PHI(:,1),neighb,xnode,model.k);

    for n = 1 : model.maxit

        PHI_aux = A*F + (I-A*K)*PHI_n;

        err = norm(PHI_aux-PHI_n,2)/norm(PHI_aux,2);

        if err >= model.tol

            PHI(:,n+1) = PHI_aux;

            PHI_n = PHI(:,n+1);

            Q=[Q,fdm2d_flux(PHI(:,n+1),neighb,xnode,model.k)];

        end

    end

    disp('Método terminado por límite de iteraciones.');

end
