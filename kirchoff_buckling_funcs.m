clear all; close all; clc

xi = sym('xi',[2,1],'real');
N1 = 0.25*(xi(1) - 1)*(xi(2) - 1); 
N2 = -0.25*(xi(1) + 1)*(xi(2) - 1);
N3 = 0.25*(xi(1) +1)*(xi(2) + 1); 
N4 = -0.25*(xi(1) - 1)*(xi(2) + 1);
%define N-matrix for the element
Ne = [N1 0 N2 0 N3 0 N4 0;
0 N1 0 N2 0 N3 0 N4];
Ne_w = [N1 N2 N3 N4];
%differentiate shape functions wrt isoparam. coordinates
dN1_dxi=gradient(N1,xi);
dN2_dxi=gradient(N2,xi);
dN3_dxi=gradient(N3,xi);
dN4_dxi=gradient(N4,xi);
%introduce node positions
xe1 = sym('xe1',[2,1],'real');
xe2 = sym('xe2',[2,1],'real');
xe3 = sym('xe3',[2,1],'real');
xe4 = sym('xe4',[2,1],'real');
%introduce spatial coordinate as fcn of isoparam. coord.
x=N1*xe1+N2*xe2+N3*xe3+N4*xe4;
%compute Jacobian
Fisop=jacobian(x,xi);
invFisop = simplify(inv(Fisop));
detFisop = simplify(det(Fisop));
P(xi) = [1, xi(1), xi(2), xi(1)^2, xi(1)*xi(2), xi(2)^2, xi(1)^3, xi(1)^2*xi(2), xi(1)*xi(2)^2, xi(2)^3, xi(1)^3*xi(2), xi(1)*xi(2)^3];
alpha = sym('alpha',[12,1],'real');
w(xi) = P(xi(1), xi(2))*alpha;
dPdxi = sym('dPdxi',[12,2],'real');
dPdxi1(xi) = diff(P(xi(1), xi(2)),xi(1));
dPdxi2(xi) = diff(P(xi(1), xi(2)),xi(2));
dPdx1(xi) = jacobian(P(xi(1),xi(2)),xi)*invFisop(:,1);
dPdx2(xi) = jacobian(P(xi(1),xi(2)),xi)*invFisop(:,2);

A = sym('A',[12,12],'real');

np = [-1 -1 ; -1 -1; -1 -1; 1 -1; 1 -1; 1 -1; 1 1; 1 1; 1 1; -1 1; -1 1; -1 1];
for i = [1 4 7 10]
    A(i,:) = P(np(i,1), np(i,2));
end
for i = [2 5 8 11]
    A(i,:) = dPdx1(np(i,1), np(i,2));
end
for i = [3 6 9 12]
    A(i,:) = dPdx2(np(i,1), np(i,2));
end

N = simplify(P*inv(A));

%matlabFunction(N, 'File','N_kirchoff_func','Vars',{xi,xe1,xe2,xe3,xe4});
dNdx = sym('dNdx',[12,2],'real');

dNdx(:,1) = diff(N(xi(1),xi(2)),xi(1))*invFisop(1,1)+diff(N(xi(1),xi(2)),xi(2))*invFisop(2,1);
dNdx(:,2) = diff(N(xi(1),xi(2)),xi(1))*invFisop(1,2)+diff(N(xi(1),xi(2)),xi(2))*invFisop(2,2);

Bast(1,:) = diff(dNdx(:,1),xi(1))*invFisop(1,1) + diff(dNdx(:,1),xi(2))*invFisop(2,1);
Bast(2,:) = diff(dNdx(:,2),xi(1))*invFisop(1,2) + diff(dNdx(:,2),xi(2))*invFisop(2,2);
Bast(3,:) = 2*(diff(dNdx(:,1),xi(1))*invFisop(1,2) + diff(dNdx(:,1),xi(2))*invFisop(2,2));

%use chain rule to compute spatial derivatives
dN1_dx=simplify( inv(Fisop)'*dN1_dxi );
dN2_dx=simplify( inv(Fisop)'*dN2_dxi );
dN3_dx=simplify( inv(Fisop)'*dN3_dxi );
dN4_dx=simplify( inv(Fisop)'*dN4_dxi );

Be = dNdx';

D = sym('D',[3,3],'real');
N_tilde = sym('N_tilde',[2,2],'real');
h = sym('h',[1,1],'real');
D_tilde = ((h^3)/12)*D;

K_ww = Bast'*D_tilde*Bast*detFisop*h;
Ge = dNdx*N_tilde*dNdx'*detFisop*h;

matlabFunction(K_ww, 'File','Kww_kirchoff_func','Vars',{xi,xe1,xe2,xe3,xe4,D,h});
matlabFunction(Ge, 'File','Ge_kirchoff_func','Vars',{xi,xe1,xe2,xe3,xe4,N_tilde,h});
