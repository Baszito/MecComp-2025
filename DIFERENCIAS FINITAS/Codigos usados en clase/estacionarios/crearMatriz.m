function [A,b,x0] = crearMatriz(N)
% TPractico Capitulos 02 y 03
% Creamos la matriz y el vector del lado derecho
% esta es la forma habitual (no sparse)

% generamos una matriz full
if 0
A=(2*diag(ones(1,N),0)-1*diag(ones(1,N-1),1)-1*diag(ones(1,N-1),-1));
A(1,[1:2])=[1 0];
A(N,[N-1:N])=[0 1];
end

% generamos una matriz sparse
if 1
unos=ones(N,1);
diagonales = [-1*unos 2*unos -1*unos];
A = spdiags(diagonales, [-1 0 1], N,N);
A(1,[1:2])=[1 0];
A(N,[N-1:N])=[0 1];
end

b = ones(N,1);
b(1)=0;
b(N)=0;


b= [b(1); (1/N^2).*ones(N-2,1); b(N)];


x0 = zeros(N,1);

endfunction