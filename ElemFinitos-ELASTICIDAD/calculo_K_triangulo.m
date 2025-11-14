#Funcion para calcular matrices K, es decir, matriz de rigidez de elementos TRIANGULARES
function [K]=calculo_K_triangulo(xnode,D)
  #Acordate en caso de tension plana multiplciar por el espesor
  #Area del triangulo, lo usamos en todas las B
  A=0.5*abs(xnode(1,1)*(xnode(2,2)-xnode(3,2)) + xnode(2,1)*(xnode(3,2)-xnode(1,2)) + xnode(3,1)*(xnode(1,2)-xnode(2,2)));

  #PRIMER NODO
  i=1; #en octave los indices arrancan en 1
  j=2;
  k=3;

  b1=xnode(j,2)-xnode(k,2);
  c1=xnode(k,1)-xnode(j,1);
  B1=zeros(3,2);
  B1=[b1 0;0 c1;c1 b1]

  #SEGUNDO NODO
  i=2;
  j=3;
  k=1;

  b2=xnode(j,2)-xnode(k,2);
  c2=xnode(k,1)-xnode(j,1);
  B2=zeros(3,2);
  B2=[b2 0;0 c2;c2 b2]

  #TERCER NODO
  i=3;
  j=1;
  k=2;

  b3=xnode(j,2)-xnode(k,2);
  c3=xnode(k,1)-xnode(j,1);
  B3=zeros(3,2);
  B3=[b3 0;0 c3;c3 b3]


  B=[B1 B2 B3];

  K=(1/(4*A)).*(B' * D * B);
  #Acordate que es : Cada una de 2x2
  #[K11 K12 K13]
  #[K21 K22 K23]
  #[K31 K32 K33]
endfunction
