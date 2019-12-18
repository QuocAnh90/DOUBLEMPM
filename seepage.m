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

%% Please select the versions of MPM!!!!!!!!!!!!!!!!!
% Original MPM: 'double_MPM'

% Please remember to change the output video file!!
Version = 'double_MPM';

%% Constutitive model
% Linear_Elastic
% Neo_Hookean_Elastic
% Mohr_Coulomb
% Water
CModel = 'Neo_Hookean_Elastic';
ppp = 0000;

%% Material porperties
E                       = 10e6           ;                  % Young modulus of solid
psp                     = 2143.0         ;                  % solid density
nu                      = 0.0            ;                  % Poison ratio
g                       = 10.0            ;                  % gravity acceleration
k                       = 1          ;
n_o                     = 0.3            ;
Kw                      = 220e6         ;
pwp                     = 1000           ;
CModel_parameter =[E,nu];

% phi                     = 40;
% psi                     = 0;
% c                       = 0000;
% CModel_parameter =[E,nu,phi,psi,c];

% K                       = 150000;                          % Bulk modulus
% u                       = 0.5;                            % vicosity
% gamma                   = 7.0;                            % gamma
% 
% CModel_parameter = [K,u,gamma];

%% Structured Grid input
NN(1)                 = 70;                              % number of nodes in X direction
NN(2)                 = 18;                             % number of nodes in Y direction
le(1)                 = 0.1;                           % size of element in X direction
le(2)                 = 0.1;                          % size of element in Y direction

%% Grid generation
[LOC,LOCC,cellCount,nodeCount] = Grid_Generation(NN,le);
% nodeCount: total number of nodes
% cellCount: total number of elements
% LOC(n,1:2): localtion coordinate of node "n" in x(1) and y(2) direction
% LOCC(e,1:2): localtion coordinate of element centroid "e" in x(1) and y(2) direction

%% Time
c                       = sqrt(Kw/pwp);
ftime                   = 10.0;
% dt                      = 0.5*le(1)/c;
dt                      = 0.00001;
ndt                     = round(ftime/dt) +1;
t                       = 0;

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

%% Particle generation
% Solid particles
particle_per_cell       = 4;
spCount                 = 280*particle_per_cell;
lsp(1)                  = le(1)/2;                                 % size of particle in X direction
lsp(2)                  = le(2)/2;                              % size of particle in Y direction
x_sp                    = zeros(spCount,1);
d_sp                    = zeros(spCount,2);

sp=1;
while sp<spCount+0.0001
    for i=1:28
        for j=1:40
            x_sp(sp,1:2)= [22*le(1)+0.5*lsp(1)+(j-1)*lsp(1) 2*le(2)+0.5*lsp(2)+(i-1)*lsp(2)];
            sp=sp+1;
        end
    end
end

% Liquid particles
wpCount                 = 240*particle_per_cell;
lwp(1)                 = le(1)/2;                                 % size of particle in X direction
lwp(2)                 = le(2)/2;                              % size of particle in Y direction
x_wp                    = zeros(wpCount,1);
d_wp                    = zeros(wpCount,2);

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
initial_figure1 = Plot_Initial(x_sp,LOC,le);
% hold on
% initial_figure2 = Plot_Initial(x_wp,LOC,le);

%% Particle variables
x_spo                   = x_sp;                                 % initial position
p_sp                    = psp * ones(spCount,1);                % Density
d_sp                    = zeros(spCount,2);                     % displacement
b_sp                    = [zeros(spCount,1) -g*ones(spCount,1)];% body force
s_sp                    = zeros(spCount,3);                     % Stress tensor
ds_sp                   = zeros(spCount,3);                     % Stress increment
v_ssp                   = zeros(spCount,2);                     % velocty
e_sp                    = zeros(spCount,3);                     % Strain tensor
de_sp                   = zeros(spCount,3);                     % Strain increment
ptraction_sp            = zeros(spCount,2);                     % traction
F_sp                    = cell(spCount,1);                      % Gradient deformation
r1_sp                   = zeros(spCount,2);
r2_sp                   = zeros(spCount,2);
n_sp                    = n_o * ones(spCount,1);
k_sp                    = k * ones(spCount,1);

x_wpo                   = x_wp;                                 % initial position
p_wp                    = psp * ones(wpCount,1);                % Density
d_wp                    = zeros(wpCount,2);                     % displacement
b_wp                    = [zeros(wpCount,1) -g*ones(wpCount,1)];% body force
pore_wp                 = zeros(wpCount,4);                     % Stress tensor
ds_wp                   = zeros(wpCount,3);                     % Stress increment
v_wp                    = zeros(wpCount,2);                     % velocty
e_wp                    = zeros(wpCount,3);                     % Strain tensor
de_wp                   = zeros(wpCount,3);                     % Strain increment
de_wwp                  = zeros(wpCount,3);                     % Strain increment
de_swp                  = zeros(spCount,3);                     % Strain increment
ptraction_wp            = zeros(wpCount,2);                     % traction
F_wp                    = cell(wpCount,1);                      % Gradient deformation
r1_wp                   = zeros(wpCount,2);
r2_wp                   = zeros(wpCount,2);

%% Initial condition
% Gradient deformation
for sp = 1:spCount
    r1_sp(sp,:) = [lsp(1,1)/2 0];
    r2_sp(sp,:) = [0 lsp(1,2)/2];
    F_sp{sp} = [1 0; 0 1];
end
r10_sp = r1_sp;
r20_sp = r2_sp;
V_sp                    = zeros(spCount,1);
for sp=1:spCount
V_sp(sp)                = 4*abs(r1_sp(sp,1)*r2_sp(sp,2)-r1_sp(sp,2)*r2_sp(sp,1)); 
end
V_spo                   = V_sp;
m_sp                    = psp * (1-n_o) * V_sp;                           % mass

for wp = 1:wpCount
    r1_wp(wp,:) = [lwp(1,1)/2 0];
    r2_wp(wp,:) = [0 lwp(1,2)/2];
    F_wp{wp} = [1 0; 0 1];
end
r10_wp = r1_wp;
r20_wp = r2_wp;
V_wp                    = zeros(wpCount,1);
for wp=1:wpCount
V_wp(wp)                = 4*abs(r1_wp(wp,1)*r2_wp(wp,2)-r1_wp(wp,2)*r2_wp(wp,1)); 
end
V_wpo                   = V_wp;
m_wp                    = pwp * n_o * V_wp;                           % mass

% Traction
for sp=1:spCount
    if x_sp(sp,2)>1.02
        ptraction_sp(sp,2) = -ppp;
    end
end

%% start the algorithm
% video
timestep = 1000;     % number of frame to save
r=timestep/20;      % number of frame per second video ~200s

writerObj2           = VideoWriter('saturatedcolumn.avi');
writerObj2.FrameRate = r;    % number of frame per second
open(writerObj2);

    for tt = 1:timestep
    ft              = ftime/timestep*tt;
%     ft=ftime;
 while t<ft+0.0000000001      
     t
     
     switch Version
         case 'double_MPM'
%         [v_ssp,x_sp,F_sp,V_sp,s_sp,p_sp,n_sp,k_sp,v_wp,x_wp,pore_wp] = double_MPM_solver(CModel,CModel_parameter,nodeCount,wpCount,spCount,cellCount,...
%         x_wp,x_sp,x_wpo,x_spo,d_wp,d_sp,le,NN,LOC,LOCC,...
%         b_sp,b_wp,V_sp,ptraction_sp,ptraction_wp,v_ssp,v_wp,s_sp,pore_wp,m_wp,m_sp,p_sp,k_sp,n_sp,...
%         nfbcx,nfbcy,fbcx,fbcy,Kw,k,F_sp,V_spo,V_wpo,n_o,dt);
    
        [v_ssp,x_sp,F_sp,V_sp,s_sp,p_sp,v_wp,x_wp,pore_wp,F_wp,V_wp] = double_MPM_solver_new2(CModel,CModel_parameter,nodeCount,wpCount,spCount,cellCount,...
        x_wp,x_sp,x_wpo,x_spo,d_wp,d_sp,F_wp,le,NN,LOC,LOCC,...
        b_sp,b_wp,V_sp,ptraction_sp,ptraction_wp,v_ssp,v_wp,s_sp,pore_wp,m_wp,m_sp,p_sp,...
        nfbcx,nfbcy,fbcx,fbcy,k,n_o,F_sp,V_spo,V_wp,dt,psp);   
    
%     [v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp] = MPM_solver(CModel,CModel_parameter,...
%     nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
%     nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt);

     end
     
     % Update time and step
        t = t+dt;
 end

 %% Plot the result
 
%   StressProfile1=Plot_Final(x_sp,LOC,le,d_sp,s_sp,v_ssp,spCount,r1_sp,r2_sp);
 
% velocity = zeros(spCount,1);
% stress = zeros(spCount,1);
% for sp=1:spCount
%     velocity(sp) = sqrt(v_ssp(sp,1)^2+v_ssp(sp,2)^2);
%     stress(sp) = s_sp(sp,2);
% end
% StressProfile1=figure;
% sz = 5;
% color = velocity;
% scatter(x_sp(:,1),x_sp(:,2),sz,color,'filled');
% grid on
% axis([0,max(LOC(:,1)),0,max(LOC(:,2))]);
% set(gca,'xtick',[0:le(1):max(LOC(:,1))]);
% set(gca,'ytick',[0:le(2):max(LOC(:,2))]);
% set(gca,'Yticklabel',[]);
% set(gca,'Xticklabel',[]);
% h=colorbar;
% colormap(jet(256))

velocity = zeros(wpCount,1);
stress = zeros(wpCount,1);
for wp=1:wpCount
    velocity(wp) = sqrt(v_wp(wp,1)^2+v_wp(wp,2)^2);
    stress(wp) = pore_wp(wp,1);
end
StressProfile1=figure('Renderer', 'painters', 'Position', [10 10 900 200]);
sz = 20;
color = velocity;
scatter(x_wp(:,1),x_wp(:,2),sz,color,'filled');
hold on
scatter(x_sp(:,1),x_sp(:,2),'g');
grid on
axis([0,max(LOC(:,1)),0,max(LOC(:,2))]);
set(gca,'xtick',[0:le(1):max(LOC(:,1))]);
set(gca,'ytick',[0:le(2):max(LOC(:,2))]);
set(gca,'Yticklabel',[]);
set(gca,'Xticklabel',[]);
h=colorbar;
colormap(jet(256))
caxis([0 5])

    frame2 = getframe(StressProfile1);
    writeVideo(writerObj2,frame2);
    end
    toc
    close(writerObj2);
