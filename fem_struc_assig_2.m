clc; clear all; close all;
xmin = 0;
xmax = 10;
ymin = 0;
ymax = 2.5;
nelx = 200;
nely = 15;
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
node_C = [];
ndofs = nnodes*5;
for node = 1:nnodes
    if coord(node,1) == xmin || coord(node,2) == ymin || coord(node,1) == xmax
        node_C = horzcat(node_C, node);
    end
end
dof_F = [1:ndofs];
dof_C = [];
for i = 1:length(node_C)
    dof_C_node = [node_C(i)*5-4 node_C(i)*5-3 node_C(i)*5-2 node_C(i)*5-1 node_C(i)*5]; % wrong order
    dof_C = horzcat(dof_C,dof_C_node);
end
a_C = zeros(1,length(dof_C))';
dof_F(dof_C) = [];

K = spalloc(ndofs,ndofs,20*ndofs); % defines K as a sparse matrix and sets the size
f = zeros(ndofs,1);
xi_v = [-1/sqrt(3) -1/sqrt(3) 1/sqrt(3) 1/sqrt(3);
        -1/sqrt(3) 1/sqrt(3) -1/sqrt(3) 1/sqrt(3)];
xi_2 = [0,0]';
pres_v = zeros(1,20);
pres_v(3:5:end) = 1;
Ex = zeros(nel,4);
Ey = zeros(nel,4);
for el = 1:nel
    fe = zeros(1,20);
    nodes = mesh(:,el);
    Ex(el,:) = coord(nodes,1);
    Ey(el,:) = coord(nodes,2); 
    y_middle = sum(coord(nodes,2))/4; 
    P = rho*g*(ymax - y_middle);
    Ke_1 = zeros(20,20);
    for i = 1:4
        xi = xi_v(:,i);
        Ke_1 = Ke_1 + Ke_mindlin_func_1(xi,coord(nodes(1),:)',coord(nodes(2),:)',coord(nodes(3),:)',coord(nodes(4),:)',D,G,h);
    end
    [Ke_2,detFisop] = Ke_mindlin_func_2(xi_2,coord(nodes(1),:)',coord(nodes(2),:)',coord(nodes(3),:)',coord(nodes(4),:)',D,G,h);
    Ke = Ke_1 + Ke_2*2;
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
                    top_trac(i*5-3) = del_t/2;
                    top_trac(j*5-3) = del_t/2;
                    fe = fe + top_trac;
                    break;
                end
            end
        end
    end
    fe = fe + detFisop*P/4*pres_v;
    K(Edof(el,2:end),Edof(el,2:end)) = K(Edof(el,2:end),Edof(el,2:end)) + Ke;
    f(Edof(el,2:end))= f(Edof(el,2:end)) + fe';
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
sfac = 1e3; % magnification factor
eldisp2(Ex,Ey,Ed_xy,plotpar,sfac);

z = 0;
for el = 1:nel
    nodes = mesh(:,el);
    [sigma(el,:), tau(el,:)] = Stress_mindlin_func_2(xi_2,coord(nodes(1),:)',coord(nodes(2),:)',coord(nodes(3),:)',coord(nodes(4),:)',Ed(el,:)',D,G,z);
    S = [sigma(el,1:2), 0, tau(el,:), sigma(el,3)]' - sum(sigma(el,1:2))/3*[1,1,1,0,0,0]';
    sigma_vm(el) = sqrt(3/2*(S'*S + sigma(el,3)^2 + tau(el,1)^2 + tau(el,2)^2));
end
figure
plotpar=[1 1 0];
sfac = 1e2; % magnification factor
eldisp2_fill(Ex,Ey,Ed_xy,plotpar,sfac,sigma_vm);
hold on
xlabel ("[m]")
ylabel ("[m]")
title("Deformed geometry for full top side displacement")
colorbar

% xi_3x3 = []
% 
% for el = 1:nel
%     nodes = mesh(:,el);
%     Ex(el,:) = coord(nodes,1);
%     Ey(el,:) = coord(nodes,2); 
%     Ke_ww = zeros(12,12);
%     for i = 1:4
%         xi = xi_v(:,i);
%         Ke_ww = Ke_ww + Kww_kirchoff_func(xi,coord(nodes(1),:)',coord(nodes(2),:)',coord(nodes(3),:)',coord(nodes(4),:)',D,G,h);
%     end
%     for i = 1:9
%         xi = xi_3x3(:,i);
%         Ge = Ge + Ge_kirchoff_func(xi,coord(nodes(1),:)',coord(nodes(2),:)',coord(nodes(3),:)',coord(nodes(4),:)',D,G,h);
%     end
%     
%     fe = fe + detFisop*P/4*pres_v;
%     K(Edof(el,2:end),Edof(el,2:end)) = K(Edof(el,2:end),Edof(el,2:end)) + Ke;
%     f(Edof(el,2:end))= f(Edof(el,2:end)) + fe';
% end



