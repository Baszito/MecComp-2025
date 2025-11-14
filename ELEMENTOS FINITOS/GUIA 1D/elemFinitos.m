function [T] = elemFinitos(xnode, model, cb, et,tip_func)
% elemFinitos: Método de Elementos Finitos 1D para conducción-reacción de calor
%
% Ecuación general:
%   rho*cp*dT/dt = k*d2T/dx2 - c*T + G
%
% Entradas:
%   xnode : vector de coordenadas nodales
%   model : struct con campos k, c, rho, cp y G
%   cb    : matriz 2x3 con condiciones de borde [tipo valor valor_extra]
%   et    : esquema temporal
%           0 = estacionario
%           1 = Euler implícito
%           2 = Euler explícito (dt crítico)
%
% Salida:
%   T     : vector de temperatura nodal final

    %--------------------------------------------------
    % 1. Parámetros del modelo
    %--------------------------------------------------
    k   = model.k;
    c   = model.c;
    rho = model.rho;
    cp  = model.cp;
    G   = model.G;

    %--------------------------------------------------
    % 2. Mallado
    %--------------------------------------------------
    Nnod = length(xnode);
    Nelem = Nnod - 1;

    %--------------------------------------------------
    % 3. Inicialización de matrices globales
    %--------------------------------------------------
    Kg = zeros(Nnod);
    Cg = zeros(Nnod);
    fg = zeros(Nnod,1);

    %--------------------------------------------------
    % 4. Ensamblado elemento a elemento
    %--------------------------------------------------
    if(tip_func==1)
      for e = 1:Nelem
          h = xnode(e+1) - xnode(e);

          % Derivadas de funciones de forma
          dN_dx = [-1, 1];

          % Matrices locales
          Ke  = (k/h) * (dN_dx' * dN_dx);
          Ce  = (h/6) * [2 1; 1 2];

          % centro del elemento
          xc = xnode(e) + h/2;
          fe  = G(xc) * (h/2) * [1; 1];

          % Ensamblado global
          idx = [e, e+1];
          Kg(idx,idx) = Kg(idx,idx) + Ke;
          Cg(idx,idx) = Cg(idx,idx) + Ce;
          fg(idx) = fg(idx) + fe;
      end
    else
      [Kg,Cg,fg,xnode]=FuncCuad(xnode,k,G);
    end
    %--------------------------------------------------
    % 5. Aplicación de condiciones de borde
    %--------------------------------------------------
    K = Kg+c*Cg;
    % --- Izquierda ---
    switch cb(1,1)
        case 1 % Dirichlet
            K(1,:) = 0; K(1,1) = 1;
            fg(1) = cb(1,2);
        case 2 % Neumann (flujo)
            fg(1) = fg(1) - cb(1,2);
        case 3 % Robin (convección)
            K(1,1) = K(1,1) + cb(1,2);
            fg(1) = fg(1) + cb(1,2) * cb(1,3);
    end

    % --- Derecha ---
    switch cb(2,1)
        case 1 % Dirichlet
            K(Nnod,:) = 0; K(Nnod,Nnod) = 1;
            fg(Nnod) = cb(2,2);
        case 2 % Neumann
            fg(Nnod) = fg(Nnod) - cb(2,2);
        case 3 % Robin
            K(Nnod,Nnod) = K(Nnod,Nnod) + cb(2,2);
            fg(Nnod) = fg(Nnod) + cb(2,2) * cb(2,3);
    end

    %--------------------------------------------------
    % 6. Resolución según esquema temporal
    %--------------------------------------------------
    % et: esquema temporal [tipo val1 val2 val3]
    % tipo==0 -> Estacionario
    % tipo==1 -> Implicito
    % tipo==2 -> Explicito
    % Si tipo==1 o tipo==2: val1=maxIt, val2=tolError
    % Si tipo==1 val3=dt
    switch et(1)
        case 0 % ESTADO ESTACIONARIO
            T = K \ fg;
            plot(xnode, T, 'b-x')

        case 1 % EULER IMPLÍCITO
            dx = min(diff(xnode));
            dt = et(4);
            Nt = et(2);
            Tn = zeros(Nnod,1); % condición inicial
            alpha = (rho * cp) / dt;
            A = inv(alpha*Cg+K);
            T_tot = Tn;
            for n = 1:Nt
                b = fg + alpha*Cg*Tn;
                T = A*b;
                T_tot = [T_tot T];

                err = (norm(T-Tn)/norm(Tn));
                plot(xnode, T, '*-b')
                title(sprintf('nit: %d - error: %e',n,err));
                pause(0.05)

                if err < et(3)
                  break
                endif
                Tn = T;
            end

        case 2
            % ------------------------------
            % EULER EXPLÍCITO
            % ------------------------------
            alpha = k / (rho * cp);
            dx = min(diff(xnode));
            dt_crit = (0.5*dx^2)/(alpha);
            dt = 0.3 * dt_crit;
            fprintf('>> dt_crit = %.4e s (usando dt=%.4e)\n', dt_crit, dt);
            Nt = et(2);
            Tn = zeros(Nnod,1);
            T_tot = Tn;
            aux = (rho * cp) / dt;
            A = aux*Cg-K;
            B = inv(Cg) * (1/aux);
            for n = 1:Nt
                T = B*(fg+A*Tn);

                # en cada iteracion corregimos/fijamos la/las
                # condiciones dirichlet
                if cb(1,1) == 1
                  T(1) = cb(1,2);
                endif
                if cb(2,1) == 1
                  T(Nnod) = cb(2,2);
                endif

                T_tot = [T_tot T];

                err = (norm(T-Tn)/norm(Tn));
                plot(xnode, T, '*-b')
                title(sprintf('nit: %d - error: %e',n,err));
                pause(0.05)

                if err < et(3)
                  break
                endif

                Tn = T;
            end

        otherwise
            fprintf('Valor de et no reconocido. Use 0 (estacionario), 1 (implícito) o 2 (explícito). \n');
            fprintf('Se resuelve el estado estacionario')
            T = K \ fg;
    end
end

