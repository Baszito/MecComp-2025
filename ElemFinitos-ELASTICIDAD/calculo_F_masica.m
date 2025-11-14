#calculo de F masica, elemento triangular}
function [Fb]=calculo_F_masica(xnode,b)
   A=0.5*abs(xnode(1,1)*(xnode(2,2)-xnode(3,2)) + xnode(2,1)*(xnode(3,2)-xnode(1,2)) + xnode(3,1)*(xnode(1,2)-xnode(2,2)));
   Fb=(A/3).*[b(1);b(2);b(1);b(2);b(1);b(2)];
endfunction
