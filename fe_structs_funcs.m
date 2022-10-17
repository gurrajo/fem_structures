clear all; close all; clc
xi = sym('xi',[2,1],'real');
N1=1-xi(1)-xi(2); N2=xi(1); N3=xi(2);
%differentiate shape functions wrt isoparam. coordinates
dN1_dxi=gradient(N1,xi);
dN2_dxi=gradient(N2,xi);
dN3_dxi=gradient(N3,xi);
%introduce node positions
xe1 = sym('xe1',[2,1],'real');
xe2 = sym('xe2',[2,1],'real');
xe3 = sym('xe3',[2,1],'real');
%introduce spatial coordinate as fcn of isoparam. coord.
x=N1*xe1+N2*xe2+N3*xe3;
%compute Jacobian
Fisop=jacobian(x,xi);
%use chain rule to compute spatial derivatives
dN1_dx=simplify( inv(Fisop)'*dN1_dxi );
dN2_dx=simplify( inv(Fisop)'*dN2_dxi );
dN3_dx=simplify( inv(Fisop)'*dN3_dxi );
%define element B-matrix
Be=[dN1_dx(1) 0 dN2_dx(1) 0 dN3_dx(1) 0;
    0 dN1_dx(2) 0 dN2_dx(2) 0 dN3_dx(2);
    dN1_dx(2) dN1_dx(1) dN2_dx(2) dN2_dx(1) dN3_dx(2) dN3_dx(1)];
matlabFunction(Be,'File','Be_cst_func','Vars',{xe1,xe2,xe3});
matlabFunction(0.5*det(Fisop),'File','Ae_cst_func','Vars',{xe1,xe2,xe3});

%define isoparametric coordinates and shape functions:
xi = sym('xi',[2,1],'real');
N1=1-xi(1)-xi(2); N2=xi(1); N3=xi(2);
%differentiate shape functions wrt isoparam. coordinates
dN1_dxi=gradient(N1,xi);
dN2_dxi=gradient(N2,xi);
dN3_dxi=gradient(N3,xi);
%introduce node positions;
Xe1 = sym('Xe1',[2,1],'real');
Xe2 = sym('Xe2',[2,1],'real');
Xe3 = sym('Xe3',[2,1],'real');
%introduce Lagrangian coord. as fcn of isoparam. coord.
X=N1*Xe1+N2*Xe2+N3*Xe3;
%compute Jacobian
Fisotr=jacobian(X,xi);
%use chain rule to compute Lagr. derivatives
dN1_dX=simplify( inv(Fisotr)'*dN1_dxi );
dN2_dX=simplify( inv(Fisotr)'*dN2_dxi );
dN3_dX=simplify( inv(Fisotr)'*dN3_dxi );
%define element B0-matrix
Be0=[dN1_dX(1) 0 dN2_dX(1) 0 dN3_dX(1) 0;
     0 dN1_dX(2) 0 dN2_dX(2) 0 dN3_dX(2);
     dN1_dX(2) 0 dN2_dX(2) 0 dN3_dX(2) 0;
     0 dN1_dX(1) 0 dN2_dX(1) 0 dN3_dX(1)];
matlabFunction(Be0,'File','Be0_cst_largedef_func','Vars',{Xe1,Xe2,Xe3});




syms Gmod lambda c2 c3 Kb 
Fv = sym('Fv',[4,1],'real');
F=[Fv(1) Fv(3) 0; Fv(4) Fv(2) 0; 0 0 1];
C=F'*F;
J=det(F);
S=Gmod*(eye(3)-inv(C))+lambda*log(J)*inv(C) + (4*c2*(trace(C)-3)+6*c3*(trace(C)-3)^2)*eye(3);
sigma=(1/J)*F*S*F';
P = F*S;
Pv=[P(1,1) P(2,2) P(1,2) P(2,1)]';
dPvdFv=sym('dPvFv',[4,4],'real');
for i=1:4
    dPvdFv(i,:)=gradient(Pv(i),Fv);
end
matlabFunction(sigma, Pv, dPvdFv,'File','yeoh_cst_func','Vars',{Gmod,lambda,c2,c3,Fv});