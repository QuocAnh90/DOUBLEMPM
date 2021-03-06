close all;
clear all;
tic;
addpath('SubFunctions/Cores');
addpath('SubFunctions/Constitutive_models');
addpath('SubFunctions/double_point');
addpath('SubFunctions/CPDI_solvers');
addpath('SubFunctions/MPM_solvers');

% Unit
% Newton - seconds - metre

%% Select algorithm here
Version = 'double_MPM';

%% Select Constutitive model here
% Linear_Elastic
% Neo_Hookean_Elastic
% Mohr_Coulomb
% Water
CModel = 'Neo_Hookean_Elastic';

%% Material porperties
E                       = 10e6           ;                  % Young modulus of solid
psp                     = 2143.0         ;                  % solid density
nu                      = 0.0            ;                  % Poison ratio
g                       = 10.0           ;                  % gravity acceleration
k                       = 1              ;                  % Permeability
n_o                     = 0.3            ;                  % Initial porosity
Kw                      = 220e6          ;                  % liquid bulk modulus
pwp                     = 1000           ;                  % Liquid density

%% Switch constitutive model
if strcmp(CModel,'Neo_Hookean_Elastic')
CModel_parameter =[E,nu];                                   % Constitutive model parameters

elseif strcmp(CModel,'Mohr_Coulomb')
phi                     = 40;                               % Friction angle
psi                     = 0;                                % Dilation angle
c                       = 0;                                % Cohesion
CModel_parameter =[E,nu,phi,psi,c];     

elseif strcmp(CModel,'Water')
K                       = 150000;                           % Bulk modulus
u                       = 0.5;                              % vicosity
gamma                   = 4;                                % gamma
CModel_parameter = [K,u,gamma];   
end

%% Structured Grid input
NN(1)                 = 70;                                 % number of nodes in X direction
NN(2)                 = 18;                                 % number of nodes in Y direction
le(1)                 = 0.1;                                % size of element in X direction
le(2)                 = 0.1;                                % size of element in Y direction

%% Grid generation
[LOC,LOCC,cellCount,nodeCount] = Grid_Generation(NN,le);
% nodeCount: total number of nodes
% cellCount: total number of elements
% LOC(n,1:2): localtion coordinate of node "n" in x(1) and y(2) direction
% LOCC(e,1:2): localtion coordinate of element centroid "e" in x(1) and y(2) direction

%% Time variables
ftime                   = 30.0;                             % Total time
dt                      = 0.00001;                          % Time step
t                       = 0;                                % Intial time step

%% Boundary nodes
% Boundary coordination
x_min = 2*le(1);
x_max = 100;
y_min = 2*le(2);
y_max = 100;
[nfbcx,nfbcy,fbcx,fbcy]=Compute_Boundary_Nodes(nodeCount,LOC,x_max,x_min,y_max,y_min);
% nfbcx: number of boundary nodes in X direction
% nfbcy: number of boundary nodes in Y direction
% fbcx: index of all boundary nodes in X direction
% fbcy: index of all boundary nodes in Y direction

%% Solid Particle generation
particle_per_cell       = 4;
spCount                 = 280*particle_per_cell;                % Total number of solid particles
lsp(1)                  = le(1)/2;                              % size of particle domain in X direction
lsp(2)                  = le(2)/2;                              % size of particle domain in Y direction
x_sp                    = zeros(spCount,1);

sp=1;
while sp<spCount+0.0001
    for i=1:28
        for j=1:40
            x_sp(sp,1:2)= [22*le(1)+0.5*lsp(1)+(j-1)*lsp(1) 2*le(2)+0.5*lsp(2)+(i-1)*lsp(2)];
            sp=sp+1;
        end
    end
end

%% Liquid Particle generation
wpCount                 = 240*particle_per_cell;                % Total number of liquid particles
lwp(1)                  = le(1)/2;                              % size of particle domain in X direction
lwp(2)                  = le(2)/2;                              % size of particle domain in Y direction
x_wp                    = zeros(wpCount,1);

wp=1;
while wp<wpCount+0.0001
    for i=1:24
        for j=1:40
            x_wp(wp,1:2)= [2*le(1)+0.5*lwp(1)+(j-1)*lwp(1) 2*le(2)+0.5*lwp(2)+(i-1)*lwp(2)];
            wp=wp+1;
        end
    end
end
% x_sp: Vector, position of MPs
% spCount: total number of MPs

%% Plot initial condition
% initial_figure1 = Plot_Initial(x_sp,LOC,le);
% hold on
% initial_figure2 = Plot_Initial(x_wp,LOC,le);

%% Solid Particle variables
x_spo                   = x_sp;                                 % initial position
p_sp                    = psp * ones(spCount,1);                % Density
d_sp                    = zeros(spCount,2);                     % displacement
b_sp                    = [zeros(spCount,1) -g*ones(spCount,1)];% body force
s_sp                    = zeros(spCount,3);                     % Stress tensor
v_ssp                   = zeros(spCount,2);                     % velocty
ptraction_sp            = zeros(spCount,2);                     % traction
F_sp                    = cell(spCount,1);                      % Gradient deformation
r1_sp                   = zeros(spCount,2);                     % CPDI Vector of particle domain in x direction
r2_sp                   = zeros(spCount,2);                     % CPDI Vector of particle domain in y direction

% Initial condition
for sp = 1:spCount
    r1_sp(sp,:) = [lsp(1,1)/2 0];        
    r2_sp(sp,:) = [0 lsp(1,2)/2];
    F_sp{sp} = [1 0; 0 1];
end

V_sp                    = zeros(spCount,1);                     % Volume
for sp=1:spCount
V_sp(sp)                = 4*abs(r1_sp(sp,1)*r2_sp(sp,2)-r1_sp(sp,2)*r2_sp(sp,1)); 
end
V_spo                   = V_sp;                                 % Initial volume
m_sp                    = psp * (1-n_o) * V_sp;           

%% Liquid Particle variables
x_wpo                   = x_wp;                                 % initial position
p_wp                    = psp * ones(wpCount,1);                % Density
d_wp                    = zeros(wpCount,2);                     % displacement
b_wp                    = [zeros(wpCount,1) -g*ones(wpCount,1)];% body force
pore_wp                 = zeros(wpCount,4);                     % Stress tensor
v_wp                    = zeros(wpCount,2);                     % velocty
ptraction_wp            = zeros(wpCount,2);                     % traction
F_wp                    = cell(wpCount,1);                      % Gradient deformation
r1_wp                   = zeros(wpCount,2);                     % CPDI Vector of particle domain in x direction
r2_wp                   = zeros(wpCount,2);                     % CPDI Vector of particle domain in y direction

% Initial condition
for wp = 1:wpCount
    r1_wp(wp,:) = [lwp(1,1)/2 0];
    r2_wp(wp,:) = [0 lwp(1,2)/2];
    F_wp{wp} = [1 0; 0 1];
end

V_wp                    = zeros(wpCount,1);                     % Volume
for wp=1:wpCount
V_wp(wp)                = 4*abs(r1_wp(wp,1)*r2_wp(wp,2)-r1_wp(wp,2)*r2_wp(wp,1)); 
end
V_wpo                   = V_wp;                                 % Initial volume
m_wp                    = pwp * n_o * V_wp;                     % mass

%% start the algorithm
% video
timestep = 3000;     % number of frame to save
r=timestep/20;       % number of frame per second video ~200s

writerObj2           = VideoWriter('seepage.avi');
writerObj2.FrameRate = r;    % number of frame per second
open(writerObj2);

    for tt = 1:timestep
    ft              = ftime/timestep*tt;
%     ft=ftime;
 while t<ft+0.0000000001      
     t
     
     switch Version
         case 'double_MPM'   
        [v_ssp,x_sp,F_sp,V_sp,s_sp,p_sp,v_wp,x_wp,pore_wp,F_wp,V_wp]...
        = double_MPM_solver(CModel,CModel_parameter,nodeCount,wpCount,spCount,cellCount,...
        x_wp,x_sp,x_wpo,x_spo,d_wp,d_sp,F_wp,le,NN,LOC,...
        b_sp,b_wp,V_sp,ptraction_sp,ptraction_wp,v_ssp,v_wp,s_sp,pore_wp,m_wp,m_sp,p_sp,...
        nfbcx,nfbcy,fbcx,fbcy,k,n_o,F_sp,V_spo,V_wp,dt,psp);   
    
%     Output
%     v_ssp       : Solid Particle velocity
%     x_sp        : Solid Particle position
%     F_sp        : Solid Particle deformation gradient
%     V_sp        : Solid Particle volume
%     s_sp        : Solid Particle Stress
%     p_sp        : Solid Particle density
%     v_wp        : Liquid Particle velocity
%     x_wp        : Liquid Particle position
%     pore_wp     : Liquid Particle pore pressure
%     F_wp        : Liquid Particle deformation gradient
%     V_wp        : Liquid Particle volume 

    % Fix the solid particles
    x_sp = x_spo;
     end
     
     % Update timestep
        t = t+dt;
 end

    %% Plot the result
    velocity = zeros(wpCount,1);
    stress = zeros(wpCount,1);
    for wp=1:wpCount
        velocity(wp) = sqrt(v_wp(wp,1)^2+v_wp(wp,2)^2);
        stress(wp) = pore_wp(wp,1);
    end
    Profile=figure('Renderer', 'painters', 'Position', [10 10 900 200]);
%     set(Profile, 'visible','off');       % Turn on off the figure plot
    sz = 20;                             % Size of particle
    color = velocity;
    scatter(x_wp(:,1),x_wp(:,2),sz,color,'filled');
    hold on
    scatter(x_sp(:,1),x_sp(:,2),2,'g');
    grid on
    axis([0,max(LOC(:,1)),0,max(LOC(:,2))]);
%     set(gca,'xtick',0:le(1):max(LOC(:,1)));
%     set(gca,'ytick',0:le(2):max(LOC(:,2)));
    yticks([0 0.5 1 1.5])
    xticks([0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5])
    grid on
    grid minor
    set(gca,'MinorGridLineStyle','-')
    h=colorbar;
    limits = [0 1];
    colormap(jet(256))
    caxis(limits)                  % Legend scale

    frame2 = getframe(Profile);
    writeVideo(writerObj2,frame2);
    end
    toc
    close(writerObj2);
