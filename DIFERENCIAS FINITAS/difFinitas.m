function [T]=difFinitas(xnode,model,cb,et)
  #------------------------------------------------------------------------#
                    #-----------------DATOS-----------------#
  #------------------------------------------------------------------------#
  #xnode es un vector de coordenadas nodales
  n=length(xnode); #tamaño de xnode
  h=xnode(2)-xnode(1); #tamaño del paso

  #separacion de las constantes del modelo :
  k=model(1); #difusividad
  c=model(2); #termino reactivo
  p=model(3); #densidad
  Cp=model(4); #calor especifico
  #NOTA : Esto considera solo fuentes constantes. Si tenes una fuente del tipo 100*x, modificar G
  G=model(5); #fuente

  #armado del sistema de ec diferenciales
  K=zeros(n,n); #matriz a resolver
  b=zeros(n,1); #matriz de coeficientes
  #------------------------------------------------------------------------#
  #--------------------------Condiciones de borde y armado de la matriz K--------------------------#
  #------------------------------------------------------------------------#
  #cb(1,1) = tipo de cond borde izq
  #DIRICHLET
  if(cb(1,1)==1)
    b(1)=cb(1,2);
    K(1,1)=1;
  #NEUMANN
  elseif(cb(1,1)==2)
    b(1)=-(G/k)*(h^2)+2*h*cb(1,2)/k;
    K(1,1)=-(2+c*h*h/k);
    K(1,2)=2;
  #ROBIN
  elseif(cb(1,1)==3)
    b(1)=-(G/k)*(h^2)-2*h*cb(1,2)*cb(1,3)/k;
    K(1,1)=-(2 + c*h*h/k + 2*h*cb(1,2)/k );
    K(1,2)=2;
  endif
  #cb(2,1) = tipo de cond
  if(cb(2,1)==1)
    b(end)=cb(2,2);
    K(n,n)=1;
  #NEUMANN
  elseif(cb(2,1)==2)
    b(end)=-(G/k)*(h^2)+2*h*cb(2,2)/k;
    K(n,n)=-(2+c*h*h/k);
    K(n,n-1)=2;
  #ROBIN
  elseif(cb(2,1)==3)
    b(n)=-(G/k)*(h^2)-2*h*cb(2,2)*cb(2,3)/k;
    K(n,n)=-(2 + c*h*h/k + 2*h*cb(2,2)/k );
    K(n,n-1)=2;
  endif

  b(2:n-1)=-(G./k)*(h^2);
  #FUNCIONALIDAD FIJA
  fila = [1 -(2+c*h*h/k) 1];
  for i=2:(n-1)
    K(i,i-1:i+1)=fila;
  end

  #et = esquema temporal
  #####CASO ESTACIONARIO
  if(et==0)
    T=K\b;
  #####FORWARD EULER
  else
    dt=0.25*h*h/k; #La mitad de los necesarios,para tener estabilidad
    T0=(zeros(size(xnode))); #condiciones iniciales nulas
    T=K\b;
    % le damos valor a la U en el nodo ficticio (t=0)
  endif



endfunction
