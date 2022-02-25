clc; clear all; close all;
xmin = 0;
xmax = 10;
ymin = 0;
ymax = 2.5;
nelx = 40; % x00 makes sure a node is placed on either end of traction boundary
nely = 10;
[mesh, coord, Edof] = rectMesh(xmin, xmax, ymin, ymax, nelx, nely); % el, u1, v1, w1, thetay1, thetax1,
E = 80*10^9;
v = 0.2;
sig_yield = 170*10^6;
p_weight = 100; %[kg]
l = 1;

h = 40*10^-3;
D = hooke(1,E,v);
G = E/(2*(1+v))*[1 , 0; 0 , 1];
rho = 1000;
g = 9.81;

total_t = -p_weight*g/l;


nel = length(mesh);
nnodes = length(coord);
node_F = [1:nnodes];
node_C = [];
ndofs = nnodes*5;
for node = 1:nnodes
    if coord(node,1) == xmin || coord(node,2) == ymin || coord(node,1) == xmax
        node_C = horzcat(node_C, node);
    end
end
node_F(node_C) = [];
dof_F = [1:ndofs];
dof_F_buckl = [1:3*nnodes];
dof_C = [];
dof_C_buckl = [];
for i = 1:length(node_C)
    dof_C_node = [node_C(i)*5-4 node_C(i)*5-3 node_C(i)*5-2 node_C(i)*5-1 node_C(i)*5]; 
    dof_C_node_buckl = [node_C(i)*3-2 node_C(i)*3-1 node_C(i)*3];
    dof_C_buckl = horzcat(dof_C_buckl,dof_C_node_buckl);
    dof_C = horzcat(dof_C,dof_C_node);
end
a_C = zeros(1,length(dof_C))';
dof_F(dof_C) = [];
dof_F_buckl(dof_C_buckl) = [];
K = spalloc(ndofs,ndofs,20*ndofs); % defines K as a sparse matrix and sets the size
f = zeros(ndofs,1);
xi_v = [-1/sqrt(3) -1/sqrt(3) 1/sqrt(3) 1/sqrt(3);
        -1/sqrt(3) 1/sqrt(3) -1/sqrt(3) 1/sqrt(3)];
xi_2 = [0,0]';

dofs = [1 2 9 14 13 3 4 10 16 15 5 6 11 18 17 7 8 12 20 19];

Ex = zeros(nel,4);
Ey = zeros(nel,4);
for el = 1:nel
   
    fe = zeros(1,20);
    nodes = mesh(:,el);
    Ex(el,:) = coord(nodes,1);
    Ey(el,:) = coord(nodes,2); 
    y_middle = sum(coord(nodes,2))/4; 

    %% Stiffness matrix
    Ke_1 = zeros(20,20);
    for i = 1:4
        xi = xi_v(:,i);
        Ke_1 = Ke_1 + Ke_mindlin_func_1(xi,coord(nodes(1),:)',coord(nodes(2),:)',coord(nodes(3),:)',coord(nodes(4),:)',D,G,h);
    end
    Ke_2 = Ke_mindlin_func_2(xi_2,coord(nodes(1),:)',coord(nodes(2),:)',coord(nodes(3),:)',coord(nodes(4),:)',D,G,h);
    Ke = Ke_1 + Ke_2;

    %% Applied traction
    for i = 1:3
        x = coord(nodes(i),1);
        y = coord(nodes(i),2);
        if x >= 4.5 && x <= 5.5 && y == 2.5
            for j = (i+1):4
                x = coord(nodes(j),1);
                y = coord(nodes(j),2);
                if x >= 4.5 && x <= 5.5 && y == 2.5
                    del_x = abs(coord(nodes(i),1)-coord(nodes(j),1));
                    del_t = total_t*del_x;
                    top_trac = zeros(1,20);
                    top_trac(i*5-3) = del_t/2; % applied in y dof for nodes on traction boundary
                    top_trac(j*5-3) = del_t/2;
                    fe = fe + top_trac;
                    break;
                end
            end
        end
    end
    %% Pressure
    P = -rho*g*(ymax - y_middle);
    for i = 1:4
        xi = xi_v(:,i);
        fe(3:5:end) = fe(3:5:end) + fe_press_mindlin_func(xi,coord(nodes(1),:)',coord(nodes(2),:)',coord(nodes(3),:)',coord(nodes(4),:)',P)';
    end

    %% Assemble
    
    K(Edof(el,2:end),Edof(el,2:end)) = K(Edof(el,2:end),Edof(el,2:end)) + Ke(dofs,dofs);

    f(Edof(el,2:end))= f(Edof(el,2:end)) + fe(dofs)';
end

a_F = K(dof_F, dof_F)\( f(dof_F) - K(dof_F, dof_C)*a_C );
f_C = K(dof_C, dof_F)*a_F + K(dof_C, dof_C)*a_C - f(dof_C); %reaction forces

a(dof_F,1) = a_F;
a(dof_C,1) = a_C;

figure
eldraw2(Ex,Ey)
figure
Ed = extract(Edof,a); % extract element displacements for plotting
Ed_xy = [Ed(:,1) , Ed(:,2) , Ed(:,6) , Ed(:,7) , Ed(:,11) , Ed(:,12) , Ed(:,16) , Ed(:,17)];
plotpar=[1 1 0];
sfac = 5e4; % magnification factor
eldisp2(Ex,Ey,Ed_xy,plotpar,sfac);

%% stress calculations
z = 0;
for el = 1:nel
    nodes = mesh(:,el);
    ae = a(Edof(el,2:end));
    [sigma(el,:), tau(el,:)] = Stress_mindlin_func(xi_2,coord(nodes(1),:)',coord(nodes(2),:)',coord(nodes(3),:)',coord(nodes(4),:)',ae,D,G,z);
    S = [sigma(el,1:2), 0, sigma(el,3),tau(el,:)]' - sum(sigma(el,1:2))/3*[1,1,1,0,0,0]';
    sigma_vm(el) = sqrt(3/2*(S'*S + sigma(el,3)^2 + tau(el,1)^2 + tau(el,2)^2));
end
figure
plotpar = [1 1 0];
sfac = 5e4; % magnification factor
eldisp2_fill(Ex,Ey,Ed_xy,plotpar,sfac,sigma_vm);
hold on
xlabel ("[m]")
ylabel ("[m]")
title("Von mises stress in deformed geometry")
colorbar

figure 
w = Ed(:,3:5:end);
fill(Ex',Ey',w');
hold on
xlabel ("[m]")
ylabel ("[m]")
title("z deflection")
colorbar

%% Task 2
gp = [-sqrt(3/5), 0, sqrt(3/5)];
wi = [5/9, 8/9];
xi_3x3 = [gp(1), gp(2), gp(3), gp(1), gp(2), gp(3), gp(1), gp(2), gp(3);
          gp(1), gp(1), gp(1), gp(2), gp(2), gp(2), gp(3), gp(3), gp(3)];
wi_3x3 = [wi(1), wi(2), wi(1), wi(1), wi(2), wi(1), wi(1), wi(2), wi(1);
          wi(1), wi(1), wi(1), wi(2), wi(2), wi(2), wi(1), wi(1), wi(1)];
Gr = spalloc(nnodes*3,nnodes*3,5*ndofs);
Kww = spalloc(nnodes*3,nnodes*3,5*ndofs);
for el = 1:nel
    Edof_Kww(el,1) = el;
    j = 2;
    for i = 1:4
        node = (Edof(el,(i*5-3)) + 4)/5; 
        Edof_Kww(el,j) = node*3-2;
        j = j+1;
        Edof_Kww(el,j) = node*3-1;
        j = j+1;
        Edof_Kww(el,j) = node*3;
        j = j+1;
    end
    nodes = mesh(:,el);
    Ex(el,:) = coord(nodes,1);
    Ey(el,:) = coord(nodes,2); 
    Ke_ww = zeros(12,12);
    Ge = zeros(12,12);
    N_tilde = h*[sigma(el,1) sigma(el,3);sigma(el,3) sigma(el,2)];
    for i = 1:4
        % w = 1 
        xi = xi_v(:,i);
        Ke_ww = Ke_ww + Kww_kirchoff_func(xi,coord(nodes(1),:)',coord(nodes(2),:)',coord(nodes(3),:)',coord(nodes(4),:)',D,h);
    end
    for i = 1:9
        xi = xi_3x3(:,i);
        w = wi_3x3(1,i)*wi_3x3(2,i);
        Ge = Ge + w*Ge_kirchoff_func(xi,coord(nodes(1),:)',coord(nodes(2),:)',coord(nodes(3),:)',coord(nodes(4),:)',N_tilde,h);
    end
    % should Ke_ww and Ge be reorganized like the mindlin Ke?
    Kww(Edof_Kww(el,2:end),Edof_Kww(el,2:end)) = Kww(Edof_Kww(el,2:end),Edof_Kww(el,2:end)) + Ke_ww;
    Gr(Edof_Kww(el,2:end),Edof_Kww(el,2:end)) = Gr(Edof_Kww(el,2:end),Edof_Kww(el,2:end)) + Ge;
end
n_lambda = 8;%
[V,D] = eigs(Kww(dof_F_buckl,dof_F_buckl),-Gr(dof_F_buckl,dof_F_buckl),n_lambda,'smallestabs');

d = diag(D);
[d_sort, indx] = sort(d);
indx_min = min(find(d_sort>0));
lambda_1 = d(indx(indx_min));
z_1 = V(:,indx(indx_min));
z_w = z_1(1:3:end);
z_w_full = zeros(nnodes,1);
z_w_full(node_F) = z_w;

figure 
fill(Ex',Ey',z_w_full(mesh));
hold on
xlabel ("[m]")
ylabel ("[m]")
title("z deflection")
colorbar
n = 0;
for i = 1:nely+1
    for j =1:nelx+1
        n = n + 1;
        X(j,i) = coord(n,1);
        Y(j,i) = coord(n,2);
        Z(j,i) = z_w_full(n);
    end
end

figure

surf(X,Y,Z)
