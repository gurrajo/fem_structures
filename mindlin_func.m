a = sym('a',[20,1],'real');
z = sym('z',[1,1],'real');
P = sym('P',[1,1],'real');
t = sym('t',[1,1],'real');
xi = sym('xi',[2,1],'real');
N1 = 0.25*(xi(1) - 1)*(xi(2) - 1); N2 = -0.25*(xi(1) + 1)*(xi(2) - 1);
N3 = 0.25*(xi(1) +1)*(xi(2) + 1); N4 = -0.25*(xi(1) - 1)*(xi(2) + 1);
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
detFisop = det(Fisop);
%use chain rule to compute spatial derivatives
dN1_dx=simplify( inv(Fisop)'*dN1_dxi );
dN2_dx=simplify( inv(Fisop)'*dN2_dxi );
dN3_dx=simplify( inv(Fisop)'*dN3_dxi );
dN4_dx=simplify( inv(Fisop)'*dN4_dxi );
%define B-matrix of element
Be_1=[dN1_dx(1), 0, dN2_dx(1), 0, dN3_dx(1), 0, dN4_dx(1), 0;
0, dN1_dx(2), 0, dN2_dx(2), 0, dN3_dx(2), 0, dN4_dx(2);
dN1_dx(2),dN1_dx(1),dN2_dx(2),dN2_dx(1),dN3_dx(2),dN3_dx(1),...
dN4_dx(2),dN4_dx(1)];

Be_2 = [dN1_dx(1), dN2_dx(1), dN3_dx(1), dN4_dx(1);
        dN1_dx(2), dN2_dx(2), dN3_dx(2), dN4_dx(2)];
D = sym('D',[3,3],'real');
G = sym('G',[2,2],'real');
h = sym('h',[1,1],'real');
D_0 = D*h;
G_0 = G*h;
D_tilde = ((h^3)/12)*D;
K_uu = detFisop*Be_1'*D_0*Be_1;
K_ww = detFisop*Be_2'*G_0*Be_2;
K_wtheta = -detFisop*Be_2'*G_0*Ne;
K_thetaw = K_wtheta';
K_thetaB = detFisop*Be_1'*D_tilde*Be_1;
K_thetaS = detFisop*Ne'*G_0*Ne;
Ke_1 = sym(zeros(20,20));
Ke_1(1:8,1:8) = K_uu;
Ke_1(13:20,13:20) = K_thetaB;


Ke_2 = sym(zeros(20,20));
Ke_2(9:12,9:12) = K_ww;
Ke_2(9:12,13:20) = K_wtheta;
Ke_2(13:20,9:12) = K_thetaw;
Ke_2(13:20,13:20) = K_thetaS;

a_u = a([1,2,6,7,11,12,16,17]);
a_w = a(3:5:end);
a_theta = a([5,4,10,9,15,14,20,19]);

sigma = D*(Be_1*a_u - z*Be_1*a_theta);

tau = G*(Be_2*a_w-Ne*a_theta);

fe_pres = Ne_w'*P*detFisop;


matlabFunction(fe_pres,'File','fe_press_mindlin_func_1','Vars',{xi,xe1,xe2,xe3,xe4,P});
matlabFunction(Ke_1,'File','Ke_mindlin_func_1','Vars',{xi,xe1,xe2,xe3,xe4,D,G,h});
matlabFunction(Ke_2, detFisop,'File','Ke_mindlin_func_2','Vars',{xi,xe1,xe2,xe3,xe4,D,G,h});
matlabFunction(sigma, tau,'File','Stress_mindlin_func_2','Vars',{xi,xe1,xe2,xe3,xe4,a,D,G,z});
