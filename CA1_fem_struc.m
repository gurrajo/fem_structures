clear all; close all; clc

B1 = 10*10^-3;
B2 = 8*10^-3;
B3 = B1*5/4;
H1 = 10*10^-3;
H2 = H1/3;
H3 = 8*10^-3;
R = 1*10^-3;
h = 100*10^-3;

rho = 0; % kg/m^3 set to zero to remove body forces
g = 9.82; % m/s^2


load('LargeDeformation_topology_coarse.mat')
[Ex,Ey] = coordxtr(Edof,Coord,Dof,3);
eldraw2(Ex,Ey)

%define size of stiffness and force
nel = length(Edof);
dofs_per_node = 2; % Number of dofs per node = 2 for 2d analysis
ndofs = max(max(Edof));
K = spalloc(ndofs,ndofs,20*ndofs); % defines K as a sparse matrix and sets the size

% to (ndof x ndof) with initial zero value
% No. of estimated non-zero entries is 20*ndofs
f = zeros(ndofs,1);
a = zeros(ndofs,1);
%% Material data
mpar.Emod = 20.e6; % Youngs modulus [Pa]
mpar.v = 0.45; % Poisson's ratio [-2
mpar.lambda = mpar.v*mpar.Emod/((1+mpar.v)*(1-2*mpar.v));
mpar.Gmod= mpar.Emod/(2*(1+mpar.v));
mpar.c2 = -mpar.Gmod/10;
mpar.c3 = mpar.Gmod/30;


ptype = 2; %ptype=1: plane stress
% 2: plane strain, 3:axisym, 4: 3d
D=hooke(ptype,mpar.Emod,mpar.v); % Constitutive matrix - plane stress
D(:,3) = [];
D(3,:) = [];
% define free dofs (dof_F) and prescribed dofs (dof_C) with their values
% (a_C)
dof_F=[1:ndofs]; 
dof_C = horzcat(dof_upper',dof_lower');
dof_F(dof_C) = []; %removing the prescribed dofs from dof_F
u_tau = 6*10^-3; % presribed displacement magnitude

%body force in each element, in x and y coordinate
grav_force = -ones(1,nel)*rho*g/3; % split evenly for each node
body = [zeros(1,nel)' , grav_force'];

% initialize stuff for newton prococedure
aold = a;
da = a;
ntime = 100; % need more times steps for finer mesh 30 ok for coarse, 100 ok for very fine
tend = 10;
t = linspace(0,tend,ntime+1);
uu = linspace(0,u_tau,ntime);

%% Setup and solve FE equations
% Assemble stiffness matrix and load vector
tol=1e-6;
for i=1:ntime
    %initial guess of unknown displacement field
    a(dof_F)=aold(dof_F);
    %%updated prescribed values
    a(dof_C) = horzcat(-1/sqrt(2)*uu(i)*ones(1,length(dof_upper)) , zeros(1,length(dof_lower)))'; %displacements vector
    %Newton iteration to find unknown displacements
    unbal=1e10; niter=0;
    while unbal > tol
        K=K.*0; f=f.*0; %nullify
        %loop over elements
         for el = 1:nel
            Ae=Ae_cst_func( [ Ex(el,1) Ey(el,1) ]',[ Ex(el,2) Ey(el,2) ]',[ Ex(el,3) Ey(el,3) ]');
            Be0=Be0_cst_largedef_func( [ Ex(el,1) Ey(el,1) ]',[ Ex(el,2) Ey(el,2) ]',[ Ex(el,3) Ey(el,3) ]');
            Fv = [1 1 0 0]' + Be0*a(Edof(el,2:end));
            
            [S,Pv,dPvdFv] = yeoh_cst_func(mpar.Gmod,mpar.lambda,mpar.c2,mpar.c3,Fv);
            
            Ke_int=Be0'*dPvdFv*Be0*Ae*h;
            fe_int=Be0'*Pv*Ae*h;
            % assembling
            K(Edof(el,2:end),Edof(el,2:end)) = K(Edof(el,2:end),Edof(el,2:end)) + Ke_int;
            f(Edof(el,2:end))= f(Edof(el,2:end)) + fe_int;
         end
        %K = real(K); % not needed when ntime is large enough
        %f = real(f); 
        %unbalance equation
        g_F = f(dof_F); %f_ext = 0
        unbal=norm(g_F);
        if unbal > tol
            %Newton update
            a(dof_F)=a(dof_F) - K(dof_F,dof_F)\g_F;
        end
        
        niter=niter+1; %update iter counter
        [niter unbal] %print on screen
        if niter>20
            disp('no convergence in Newton iteration')
            break;
        end
    end
    %save data for post-processing and update state variables now
    % (when Newton iteration has converged in the time increment)
    %save how large the displacement has changed during the timesstep
    % to get a good initial guess for next timestep
    f_C = K(dof_C, dof_F)*a(dof_F) + K(dof_C, dof_C)*a(dof_C) - f(dof_C); %reaction forces
    Q = zeros(size(f));
    Q(dof_C) = f_C;
    top_force_hori(i) = sum(Q(dof_upper(1:2:end)));
    top_force_vert(i) = sum(Q(dof_upper(2:2:end)));
    
    da=a-aold;
    aold = a;
    if mod(i,10) == 0
        % plot every 10th timestep
        figure
        Ed = extract(Edof,a); % extract element displacements for plotting
        plotpar=[1 1 0];
        sfac = 1; % magnification factor
        eldisp2(Ex,Ey,Ed,plotpar,sfac);
    end
end
figure
plot(uu,top_force_hori)
hold on
plot(uu,top_force_vert)
legend("horizontal force","vertical force")
xlabel("u_\tau [m]")
ylabel("force [N]")
%find the stresses in the elements sigx, sigy, tauxy

for el = 1:nel
    Be0=Be0_cst_largedef_func( [ Ex(el,1) Ey(el,1) ]',[ Ex(el,2) Ey(el,2) ]',[ Ex(el,3) Ey(el,3) ]');
    Fv = [1 1 0 0]' + Be0*a(Edof(el,2:end));
    [S,~,~] =  yeoh_cst_func(mpar.Gmod,mpar.lambda,mpar.c2,mpar.c3,Fv);
    Sdev=S-sum(diag(S))/3*eye(3);
    vonMises_stress(el)=sqrt(3/2)*sqrt(sum(sum(Sdev.*Sdev)));
end

vert_disp = a(2:2:end);
disp("maximum vertical displacemnt: " + max(abs(vert_disp))*1e3 + " [mm]")

figure
fill(Ex',Ey',vonMises_stress)
title('\sigma_v_o_n [Pa]')
colorbar
axis equal

