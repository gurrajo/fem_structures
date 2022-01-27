clear all; close all; clc
%define isoparametric coordinates and shape functions:
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

Z = 30; 
H = 15;
B = 10;
D = 4;
b = 4;
h_len = 6;
r = 1; %[m]

rho = 2300; % kg/m^3
g = 9.82; % m/s^2

rho_vat = 997; % kg/m^3
p = -Z*rho_vat*g;

Nr = 5 ; 
Nt = 20 ;
elemtype = 1;

[Edof,Coord,Ex,Ey,LeftSide_nodes,TopSide_nodes,RightSide_nodes,BottomSide_nodes] = TunnelMeshGen(H,B,D,b,h_len,r,Nr,Nt,elemtype);
eldraw2(Ex,Ey)

%define size of stiffness and force
nel = length(Edof);
dofs_per_node = 2; % Number of dofs per node = 2 for 2d analysis
ndofs = max(max(Edof));
K = spalloc(ndofs,ndofs,20*ndofs); % defines K as a sparse matrix and sets the size

Dof = [1:ndofs];
% to (ndof x ndof) with initial zero value
% No. of estimated non-zero entries is 20*ndofs
f = zeros(ndofs,1);
a = zeros(ndofs,1);
%% Material data
mpar.Emod = 40.e9; % Youngs modulus [Pa]
mpar.v = 0.1; % Poisson's ratio [-]

ptype = 2; %ptype=1: plane stress
% 2: plane strain, 3:axisym, 4: 3d
D=hooke(ptype,mpar.Emod,mpar.v); % Constitutive matrix - plane stress
D(:,3) = [];
D(3,:) = [];
% define free dofs (dof_F) and prescribed dofs (dof_C) with their values
% (a_C)
dof_F=[1:ndofs]; 
dof_C = horzcat(LeftSide_nodes'*2 - ones(1, length(LeftSide_nodes)) , BottomSide_nodes'*2, BottomSide_nodes'*2 - ones(1, length(BottomSide_nodes)));
dof_F(dof_C) = []; %removing the prescribed dofs from dof_F
a_C = zeros(1,length(LeftSide_nodes) + length(BottomSide_nodes)*2)'; %prescribed displacements
%body force in each element, in x and y coordinate
grav_force = -ones(1,nel)*rho*g/3; % split evenly for each node
body = [zeros(1,nel)' , grav_force'];
%% Setup and solve FE equations
% Assemble stiffness matrix and load vector
h = 1;
for el = 1:nel
    
    Ae=Ae_cst_func( [ Ex(el,1) Ey(el,1) ]',[ Ex(el,2) Ey(el,2) ]',[ Ex(el,3) Ey(el,3) ]');
    Be=Be_cst_func( [ Ex(el,1) Ey(el,1) ]',[ Ex(el,2) Ey(el,2) ]',[ Ex(el,3) Ey(el,3) ]');
    Ke=Be'*D*Be*Ae*h;
    fe=[body(el,:) body(el,:) body(el,:)]'*Ae*h;

    nodes = Edof(el,3:2:7)/2;
    is_top_nodes = ismember(nodes, TopSide_nodes);
    is_right_nodes = ismember(nodes, RightSide_nodes);
    if sum(is_top_nodes) > 1
        top_nodes = nodes(is_top_nodes);
        del_x = abs(diff(Ex(el,is_top_nodes)));
        force = zeros(1,6);
        is_top_dofs = zeros(1,6);
        is_top_dofs(2:2:6) = is_top_nodes;
        force(logical(is_top_dofs)) = del_x*p/2;
        fe = fe + force';
    elseif sum(is_right_nodes) > 1
        right_nodes = nodes(is_right_nodes);
        del_y = abs(diff(Ey(el,is_right_nodes)));
        is_right_dofs = zeros(1,6);
        force = zeros(1,6);
        is_right_dofs(1:2:5) = is_right_nodes;
        force(logical(is_right_dofs)) = del_y*p/2;
        fe = fe + force';
    end
    % assembling
    K(Edof(el,2:end),Edof(el,2:end)) = K(Edof(el,2:end),Edof(el,2:end)) + Ke;
    f(Edof(el,2:end))= f(Edof(el,2:end)) + fe;
end
Ks=sparse(K); %use the sparse structure of K
a_F = Ks(dof_F, dof_F)\( f(dof_F) - Ks(dof_F, dof_C)*a_C );
f_C = Ks(dof_C, dof_F)*a_F + Ks(dof_C, dof_C)*a_C - f(dof_C); %reaction forces

a(dof_F,1) = a_F;
a(dof_C,1) = a_C;
Q = zeros(size(f));
Q(dof_C) = f_C;

figure
Ed = extract(Edof,a); % extract element displacements for plotting
plotpar=[1 1 0];
sfac = 1e4; % magnification factor
eldisp2(Ex,Ey,Ed,plotpar,sfac);
%find the stresses in the elements sigx, sigy, tauxy
Es = zeros(nel,3);Et = zeros(nel,3);
for el = 1:nel
    Be=Be_cst_func( [ Ex(el,1) Ey(el,1) ]',[ Ex(el,2) Ey(el,2) ]',[ Ex(el,3) Ey(el,3) ]');

    Et(el,:)=Be*a(Edof(el,2:end));
    Es(el,:)=D*Be*a(Edof(el,2:end));
end
%for plotting of stress in element (here: Es(el,1)=sigma_xx)
for el= 1:nel
    Esm(el,1:3)=ones(1,3)*Es(el,1);
    Esm_2(el,1:3) = ones(1,3)*Es(el,2);
end

figure
fill(Ex',Ey',Esm')
title('\sigma_x_x [Pa]')
colorbar
axis equal

figure
fill(Ex',Ey',Esm_2')
title('\sigma_y_y stress [Pa]')
colorbar
axis equal

vert_disp = a(2:2:end);
disp("maximum vertical displacemnt: " + max(abs(vert_disp))*1e3 + " [mm]")

% mu and lambda
mu = mpar.Emod/(2*(1+mpar.v));
lambda = mpar.v*mpar.Emod/((1 + mpar.v)*(1-2*mpar.v));

for i = 1:nel
    sig_z  = lambda*(Es(i,1) + Es(i,2))/(2*(mu + lambda));
    S = [Es(i,1) Es(i,3) 0 ; Es(i,3) Es(i,2) 0 ; 0 0 sig_z];
    principal_stress(i,1) = max(eig(S));

end

figure
fill(Ex',Ey',principal_stress')
title('\sigma_m_a_x [Pa]')
colorbar
axis equal