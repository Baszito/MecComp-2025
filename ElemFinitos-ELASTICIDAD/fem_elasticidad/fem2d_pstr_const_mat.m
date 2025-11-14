function [D] = fem2d_pstr_const_mat(model)
% Descripción: módulo para calcular la matriz constitutiva del sistema. La misma 
% depende de los valores del Módulo de Young y el Coeficiente de Poisson. La forma
% de esta matriz es siempre de 3x3, pero la selección de coeficientes determina el
% tratamiento del sistema como Deformación o Tensión Plana

% Entrada:
% * model: struct con todos los datos del modelo (constantes, esquema numérico, etc.)

% Salida:
% * D: matriz constitutiva del sistema.
% ----------------------------------------------------------------------
    D=zeros(3,3);
    E=model.young;
    v=model.poiss;
  if(model.pstrs==1) %TENSION PLANA
    D(1,1)=E/(1-v^2);
    D(2,2)=E/(1-v^2);
    D(1,2)=v*D(1,1);
    D(2,1)=v*D(1,1);
    D(3,3)=E/(2*(1+v));
  else
    D(1,1)=(E*(1-v))/((1+v)*(1-2*v));
    D(2,2)=(E*(1-v))/((1+v)*(1-2*v));
    D(1,2)=(v/(1-v))*D(1,1);
    D(2,1)=(v/(1-v))*D(1,1);
    D(3,3)=E/(2*(1+v));
  end
end