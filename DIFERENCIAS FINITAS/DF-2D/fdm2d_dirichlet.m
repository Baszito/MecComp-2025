function [K,F] = fdm2d_dirichlet(K,F,DIR)
    M = size(DIR, 1);
    for n = 1 : M
    % 1) Rellenar la fila DIR(n) de K con 0 y en la columna central de dicha
    % fila colocar un uno.
        K(DIR(n),:)=zeros(length(K(1,:)),1);
        K(DIR(n),DIR(n))=1;
    % 2) Asignar en F(DIR(n)) el valor de la condicion Dirichlet 
        F(DIR(n,1))=DIR(n,2);
    end
end

