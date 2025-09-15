function [T]=difFinitas(xnode,model,cb,et)
  #SOLO ESQUEMAS CON MALLAS HOMOGENEAS
  #------------------------------------------------------------------------#
                    #-----------------DATOS-----------------#
  #------------------------------------------------------------------------#
  #xnode es un vector de coordenadas nodales
  n=length(xnode); #tamaño de xnode
  dx=xnode(2)-xnode(1); #tamaño del paso

  #separacion de las constantes del modelo :
  k=model.k; #difusividad
  c=model.c; #termino reactivo
  rho=model.rho; #densidad
  cp=model.cp; #calor especifico
  #NOTA : Esto considera solo fuentes constantes. Si tenes una fuente del tipo 100*x, modificar Q
  Q=model.Q; #fuente

  #armado del sistema de ec diferenciales
  K=zeros(n,n); #matriz a resolver
  b=zeros(n,1); #matriz de coeficientes
  #FUNCIONALIDAD FIJA
  fila = [1 -(2+c*dx*dx/k) 1];
  for i=2:(n-1)
    K(i,i-1:i+1)=fila;
    b(i)=(dx^2)*Q(i)/k;
  end
  #------------------------------------------------------------------------#
  #--------------------------Condiciones de borde y armado de la matriz K--------------------------#
  #------------------------------------------------------------------------#
  #cb(1,1) = tipo de cond borde izq
  #NEUMANN
  if(cb(1,1)==2)
    b(1)=(Q(1)/k)*(dx^2) - 2*dx*cb(1,2)/k;
    K(1,1)=-(2+c*dx*dx/k);
    K(1,2)=2;
  #ROBIN
  elseif(cb(1,1)==3)
    b(1)=(Q(1)/k)*(dx^2) + 2*dx*cb(1,2)*cb(1,3)/k;
    K(1,1)=-(2 + c*dx*dx/k + 2*dx*cb(1,2)/k );
    K(1,2)=2;
  endif
  #cb(2,1) = tipo de cond
  #NEUMANN
  if(cb(2,1)==2)
    b(end)=(Q(n)/k)*(dx^2)+2*dx*cb(2,2)/k;
    K(n,n)=-(2+c*dx*dx/k);
    K(n,n-1)=2;
    #ROBIN
    elseif(cb(2,1)==3)
      b(n)=(Q(n)/k)*(dx^2) + 2*dx*cb(2,2)*cb(2,3)/k;
      K(n,n)=-(2 + c*dx*dx/k + 2*dx*cb(2,2)/k );
      K(n,n-1)=2;
endif
  #et = esquema temporal
  #####CASO ESTACIONARIO
  if(et==0)
    b*=-1;
    if(cb(1,1)==1)
      b(1)=cb(1,2);
      K(1,1)=1;
    endif
    if(cb(2,1)==1)
      b(n)=cb(2,2);
      K(n,n)=1;
    endif
    T=K\b;
    plot(xnode,T,'*-');
    grid on;
  else
  IDENTIDAD=eye(n);
  T=zeros(length(xnode),1);
  tol=0.001; #Cuando la dif entre un paso y otro es menor a TOL, se frena el for
  err=100000;
  maxit=1000;
  it=0;
  dt=0.25*dx*dx/k;
  if(et==1)#Forward euler
    while((err>tol) & (it<maxit))
    it+=1;
    figure(1)
    plot(xnode,T,'*-')
    title(sprintf('t = %5.3f',it*dt));
    grid on;
    grid minor;
    pause(0.01);
      T_aux=T;
      T=((dt/rho*cp)*(K*T_aux+b))+T_aux;
      if(cb(1,1)==1)
        T(1)=cb(1,2);#CONDICION DIRICHLET
      endif
      if(cb(2,1)==1)
        T(end)=cb(2,2);#CONDICION DIRICHLET
      endif
      err=abs(T(2:n-1)-T_aux(2:n-1));
     endwhile
    elseif(et==2)#Backward euler
    K_inv=inv(((rho*cp/dt)*IDENTIDAD-K));
    while((err>tol) & (it<maxit))
    it+=1;
    figure(1)
    plot(xnode,T,'*-')
    title(sprintf('t = %5.3f',it*dt));
    grid on;
    grid minor;
    pause(0.01);
    T_aux=T;
     T=K_inv*(b+(rho*cp/dt)*T);
      if(cb(1,1)==1)
        T(1)=cb(1,2);#CONDICION DIRICHLET
      endif
      if(cb(2,1)==1)
        T(end)=cb(2,2);#CONDICION DIRICHLET
      endif
      err=abs(T(2:n-1)-T_aux(2:n-1));
     endwhile
    elseif(et==3)#Crank-Nicholson
    K_inv=inv((rho*cp/dt)*IDENTIDAD - 0.5*K);
    while((err>tol) & (it<maxit))
    it+=1;
    figure(1)
    plot(xnode,T,'*-')
    title(sprintf('t = %5.3f',it*dt));
    grid on;
    grid minor;
    pause(0.01);
    T_aux=T;
     T=K_inv*(((0.5*K + (rho*cp/dt)*IDENTIDAD)*T) +b);
      if(cb(1,1)==1)
        T(1)=cb(1,2);#CONDICION DIRICHLET
      endif
      if(cb(2,1)==1)
        T(end)=cb(2,2);#CONDICION DIRICHLET
      endif
      err=abs(T(2:n-1)-T_aux(2:n-1));
     endwhile
   endif
  endif

endfunction
