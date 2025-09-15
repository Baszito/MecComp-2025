#EjerParcial DIF FINITAS


clc;
clear all;
xnode=-4:4/10:4;
model.k=3;
model.c=0;
model.Q=100.*xnode;
model.rho=0;
model.cp=0;
cb=[1 20 0;3 2 10];
et=0;

[T_dir]=difFinitas(xnode,model,cb,et)
