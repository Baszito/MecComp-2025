function Q=cuad_gauss_c(f,a,b,L,n)
% f: funcion a calcular
% a y b limites de integracion
% L: no de intervaloes
% n=grado
[xg,w]=gauss_xw(n);

x=linspace(a,b,L+1);
h=(b-a)/L;
Q=0;
for i=1:L
  lambdat=h/2*(xg+1)+x(i);
  Q+=h/2*(w'*f(lambdat));
endfor
