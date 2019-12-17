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
CModel = 'Linear_Elastic';
ppp = 0000;

%% Material porperties
E                       = 10e6           ;                  % Young modulus of solid
psp                     = 2143.0         ;                  % solid density
nu                      = 0.0            ;                  % Poison ratio
g                       = 10.0            ;                  % gravity acceleration
k                       = 0.001          ;
n_o                     = 0.5            ;
Kw                      = 2200e6         ;
pwp                     = 1000           ;
CModel_parameter =[E,nu];

%% Structured Grid input
NN(1)                 = 55;                              % number of nodes in X direction
NN(2)                 = 6;                             % number of nodes in Y direction
le(1)                 = 0.02;                           % size of element in X direction
le(2)                 = 0.02;                          % size of element in Y direction

%% Grid generation
[LOC,LOCC,cellCount,nodeCount] = Grid_Generation(NN,le);
% nodeCount: total number of nodes
% cellCount: total number of elements
% LOC(n,1:2): localtion coordinate of node "n" in x(1) and y(2) direction
% LOCC(e,1:2): localtion coordinate of element centroid "e" in x(1) and y(2) direction

%% Time
c                       = sqrt(E/psp);
ftime                   = 0.2;
dt                      = 0.000001;
ndt                     = round(ftime/dt) +1;
t                       = 0;

%% Boundary nodes
% Boundary coordination
x_min = 2*le(1);
x_max = 100;
y_min = 2*le(2);
y_max = 2*le(2);
[nfbcx,nfbcy,fbcx,fbcy]=Compute_Boundary_Nodes(nodeCount,LOC,x_max,x_min,y_max,y_min);
% nfbcx: number of boundary nodes in X direction
% nfbcy: number of boundary nodes in Y direction
% fbcx: index of all boundary nodes in X direction
% fbcy: index of all boundary nodes in Y direction

%% Particle generation
% Solid particles
particle_per_cell       = 4;
spCount                 = 20*particle_per_cell;
lsp(1,1)                 = le(1)/2;                                 % size of particle in X direction
lsp(1,2)                 = le(2)/2;                              % size of particle in Y direction
x_sp                    = zeros(spCount,1);
d_sp                    = zeros(spCount,2);

sp=1;
while sp<spCount+0.0001
    for i=1:2
        for j=1:40
            x_sp(sp,1:2)= [2*le(1,1)+0.5*lsp(1,1)+(j-1)*lsp(1,1) 2*le(1,2)+0.5*lsp(1,2)+(i-1)*lsp(1,2)];
            sp=sp+1;
        end
    end
end

% Liquid particles
wpCount                 = 50*particle_per_cell;
lwp(1,1)                 = le(1)/2;                                 % size of particle in X direction
lwp(1,2)                 = le(2)/2;                              % size of particle in Y direction
x_wp                    = zeros(wpCount,1);
d_wp                    = zeros(wpCount,2);

wp=1;
while wp<wpCount+0.0001
    for i=1:2
        for j=1:100
            x_wp(wp,1:2)= [2*le(1,1)+0.5*lwp(1,1)+(j-1)*lwp(1,1) 2*le(1,2)+0.5*lwp(1,2)+(i-1)*lwp(1,2)];
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
% b_sp                    = [zeros(spCount,1) -g*ones(spCount,1)];% body force
b_sp                    = [-g*ones(spCount,1) zeros(spCount,1)];% body force
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
n_spo                    = n_sp;
k_sp                    = k * ones(spCount,1);

x_wpo                   = x_wp;                                 % initial position
p_wp                    = psp * ones(wpCount,1);                % Density
d_wp                    = zeros(wpCount,2);                     % displacement
% b_wp                    = [zeros(wpCount,1) -g*ones(wpCount,1)];% body force
b_wp                    = [-g*ones(wpCount,1) zeros(wpCount,1)];% body force

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

for wp = 1:wpCount
    if x_wpo(wp,1)> 0.44
        n_wp(wp)                    = 1;
    else
        n_wp(wp)                    = n_o;
    end   
    n_wpo                    = n_wp;
end
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

m_wp(wp)                    = pwp * n_wpo(wp) * V_wp(wp);                           % mass
end

V_wpo                   = V_wp;

% Traction
for sp=1:spCount
    if x_sp(sp,2)>1.02
        ptraction_sp(sp,2) = -ppp;
    end
end

% Pore water pressure
PORE = cell(wpCount,1);
for wp=1:wpCount
    PORE{wp} = [-ppp 0; 0 -ppp];
end

%% start the algorithm
% video
timestep = 200;     % number of frame to save
r=timestep/20;      % number of frame per second video ~200s

writerObj2           = VideoWriter('2.avi');
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
    
     end
     
     % Update time and step
        t = t+dt;
 end

 %% Plot the result
     for wp=1:wpCount
        x_wp(wp,2)=x_wp(wp,2)-2*le(2);
     end
     
     for sp=1:spCount
        x_sp(sp,2)=x_sp(sp,2)-2*le(2);
     end
     
   StressProfile1=figure;
%     set(StressProfile1, 'visible','off');
%     plot(x_sp(:,1),s_sp(:,1),'.');
    plot(-s_sp(:,1),x_sp(:,1),'x',-pore_wp(:,1),x_wp(:,1),'o');
%     plot(-pore_wp,x_wp(:,2),'o');
    ylabel('Length (Pa)'); % label for y axis
    xlabel('Stress (m)'); % label for x axis
    legend('effective stress','pore pressure','Location','northeast')
%     axis([0 ppp*1.25 0 1.02]);
%     axis([0 20000 0 1]);
    
         for wp=1:wpCount
        x_wp(wp,2)=x_wp(wp,2)+2*le(2);
         end
         
         for sp=1:spCount
        x_sp(sp,2)=x_sp(sp,2)+2*le(2);
        end

     
    frame2 = getframe(StressProfile1);
    writeVideo(writerObj2,frame2);
    end
    toc
    close(writerObj2);
