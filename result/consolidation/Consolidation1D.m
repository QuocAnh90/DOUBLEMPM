% close all;
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
ppp = 10000;

%% Material porperties
E                       = 10e6           ;                  % Young modulus of solid
psp                     = 2143.0         ;                  % solid density
nu                      = 0.0            ;                  % Poison ratio
g                       = 0.0            ;                  % gravity acceleration
k                       = 0.001          ;
n_o                     = 0.3            ;
Kw                      = 2200e6         ;
pwp                     = 1000           ;
CModel_parameter =[E,nu];

%% Structured Grid input
NN(1)                 = 6;                              % number of nodes in X direction
NN(2)                 = 55;                             % number of nodes in Y direction
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
time                    = [];
error                   = [];

%% Boundary nodes
% Boundary coordination
x_min = 0.04;
x_max = 0.06;
y_min = 0.04;
y_max = 1.06;
[nfbcx,nfbcy,fbcx,fbcy]=Compute_Boundary_Nodes(nodeCount,LOC,x_max,x_min,y_max,y_min);
% nfbcx: number of boundary nodes in X direction
% nfbcy: number of boundary nodes in Y direction
% fbcx: index of all boundary nodes in X direction
% fbcy: index of all boundary nodes in Y direction

%% Particle generation
% Solid particles
particle_per_cell       = 4;
spCount                 = 50*particle_per_cell;
lsp(1,1)                 = le(1)/2;                                 % size of particle in X direction
lsp(1,2)                 = le(2)/2;                              % size of particle in Y direction
x_sp                    = zeros(spCount,1);
d_sp                    = zeros(spCount,2);

sp=1;
while sp<spCount+0.0001
    for i=1:100
        for j=1:2
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
    for i=1:100
        for j=1:2
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

x_wpo                   = x_sp;                                 % initial position
p_wp                    = psp * ones(spCount,1);                % Density
d_wp                    = zeros(spCount,2);                     % displacement
b_wp                    = [zeros(wpCount,1) -g*ones(wpCount,1)];% body force
pore_wp                 = zeros(spCount,1);                     % Stress tensor
ds_wp                   = zeros(spCount,3);                     % Stress increment
v_wp                    = zeros(spCount,2);                     % velocty
e_wp                    = zeros(spCount,3);                     % Strain tensor
de_wp                   = zeros(spCount,3);                     % Strain increment
de_wwp                  = zeros(spCount,3);                     % Strain increment
de_swp                  = zeros(spCount,3);                     % Strain increment
ptraction_wp            = zeros(spCount,2);                     % traction
F_wp                    = cell(spCount,1);                      % Gradient deformation
r1_wp                   = zeros(spCount,2);
r2_wp                   = zeros(spCount,2);

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
V_wp(wp)                = 4*abs(r1_sp(wp,1)*r2_wp(wp,2)-r1_wp(wp,2)*r2_wp(wp,1)); 
end
V_wpo                   = V_wp;
m_wp                    = pwp * n_o * V_wp;                           % mass

% Traction
for sp=1:spCount
    if x_sp(sp,2)>1.02
        ptraction_sp(sp,2) = -ppp;
    end
end

% Pore water pressure
for wp=1:wpCount
    pore_wp(wp) = -ppp;
end

PORE = cell(wpCount,1);
for wp=1:wpCount
    PORE{wp} = [pore_wp(wp) 0; 0 pore_wp(wp)];
end

%% start the algorithm
% video
timestep = 200;     % number of frame to save
r=timestep/20;      % number of frame per second video ~200s

writerObj2           = VideoWriter('consolidation.avi');
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
    
        [v_ssp,x_sp,F_sp,V_sp,s_sp,p_sp,v_wp,x_wp,pore_wp,V_wp] = double_MPM_solver_new1(CModel,CModel_parameter,nodeCount,wpCount,spCount,cellCount,...
        x_wp,x_sp,x_wpo,x_spo,d_wp,d_sp,le,NN,LOC,LOCC,...
        b_sp,b_wp,V_sp,ptraction_sp,ptraction_wp,v_ssp,v_wp,s_sp,pore_wp,m_wp,m_sp,p_sp,...
        nfbcx,nfbcy,fbcx,fbcy,Kw,k,F_sp,V_spo,V_wp,dt,psp);

    
     end
     
     % Update time and step
        t = t+dt;
        
 end

 %% Plot the result
 
  % Store Tv
       mv =	1/E;
       Sp=0.3*(1/Kw)+0.7/E;
       k_ana = mean(k_sp);
    
       Cv = k_ana/(mv+Sp)/pwp;
       Tv=Cv*t;
       L = 1;
       
 pore_wp_ana = zeros(wpCount,1);
 
     for wp=1:wpCount
        x_wp(wp,2)=x_wp(wp,2)-2*le(2);
     end
     
     for sp=1:spCount
        x_sp(sp,2)=x_sp(sp,2)-2*le(2);
     end
     
     % analytical solution of the effective and pore water pressure
 for i=1:100
     M_i = (i-0.5)*pi;
     for wp=1:wpCount
         pore_wp_ana(wp) = pore_wp_ana(wp) + 2/M_i*sin(M_i*(abs(x_wp(wp,2)-1))/(L))*exp(-M_i*M_i*Tv);
     end
 end
 
 for wp=1:wpCount
 pore_wp_ana(wp)=1/(mv+Sp)  *  (log(1+(exp((mv+Sp)*ppp)-1)*pore_wp_ana(wp)));
 end

    StressProfile1=figure;
    set(StressProfile1, 'visible','off');
%     plot(x_sp(:,1),s_sp(:,1),'.');
    plot(-s_sp(:,2),x_sp(:,2),'x',pore_wp,x_wp(:,2),'o',pore_wp_ana,x_wp(:,2));
    ylabel('Length (Pa)'); % label for y axis
    xlabel('Stress (m)'); % label for x axis
    legend('effective stress','pore pressure','analytical pore pressure','Location','northeast')
%     axis([0 ppp*1.25 0 1.02]);
    title('consolidation')
         for wp=1:wpCount
        x_wp(wp,2)=x_wp(wp,2)+2*le(2);
         end
         
         for sp=1:spCount
        x_sp(sp,2)=x_sp(sp,2)+2*le(2);
        end

     time = [time t];
     error = [error norm(pore_wp-pore_wp_ana)];
     
    frame2 = getframe(StressProfile1);
    writeVideo(writerObj2,frame2);
    end
    toc
    close(writerObj2);
