#Funcion para armar matriz D en caso de deformacion plana o tension plana
function [D]=armar_D(E,v,tipo)
  D=zeros(3,3);
  if(tipo==0) #TENSION PLANA
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
  endif
endfunction

