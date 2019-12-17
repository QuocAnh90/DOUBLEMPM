function[v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp,r1_sp,r2_sp,a_sp]...
    = CPDI_solver_alpha(CModel,CModel_parameter,...
    nodeCount,spCount,cellCount,...
    x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,...
    v_ssp,s_sp,m_sp,p_sp,nfbcx,nfbcy,fbcx,fbcy,...
    F_sp,V_spo,dt,r1_sp,r10_sp,r2_sp,r20_sp,a_sp)
     
%% Store particles into cell
[N,dN,CONNECT,spElems,mspoints,NODES] = Compute_Interpolator_CPDI(spCount,cellCount,x_sp,le,NN,LOC,r1_sp,r2_sp,V_sp);
% N: Shape function
% dN: Gradient of shape function
% spElems(p): element index where "p" locates
% mspoints(e): all particles indexes where locate in the cell "e"

 %% Mapping from particle to nodes
[nmass_si,nmomentum_si,niforce_si,neforce_si,traction_si]=Interpolate_Particle_To_Grid(NODES,nodeCount,CONNECT,le,N,dN,spCount,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp);

%% Interpolation from particle to grid task
% Node variables
nacc_old                 = zeros(nodeCount,2);

 for p=1:spCount
 for j=1:NODES(p)
     npid                           = CONNECT{p}(j);
     
          if N{p}(j)==1
         continue
          end
          
 % Momentum
 nacc_old(npid,:)      = nacc_old(npid,:) + m_sp(p)*a_sp(p,:)*N{p}(j);
 end 
 end
 
 for n = 1:nodeCount
    if nmass_si(n) ~=0
    nacc_old(n,:)       = nacc_old(n,:)/nmass_si(n);
    end
 end

 % Compute nodal velocity
nvelo_si = zeros(nodeCount,2);
for n = 1:nodeCount
    if nmass_si(n) ~=0
    nvelo_si(n,:)   = nmomentum_si(n,:)/nmass_si(n);
    end
end
[nvelo_si] = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nvelo_si); % Boundary condition for nodal force


%% Update momentum
% Update force and momentum
 nforce_si      = niforce_si + neforce_si + traction_si;
%  nmomentum_si   = nmomentum_si + nforce_si*dt;
 
 [nforce_si]  = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nforce_si); % Boundary condition for nodal force
% [nmomentum_si]  = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nmomentum_si); % Boundary condition for nodal force

 % time integration
p_b = 0.818;
a_m = (2*p_b-1)/(1+p_b);
beta = (5-3*p_b)/(1+p_b).^2/(2-p_b);
gamma = 3/2-a_m;

 nacc_middle = zeros(nodeCount,2);
 nacc = zeros(nodeCount,2);
 
 for n = 1:nodeCount
     if nmass_si(n) ~=0
    nacc_middle(n,:)    = nforce_si(n,:)/nmass_si(n);
    nacc(n,:)           = (nacc_middle(n,:) - a_m *  nacc_old(n,:))/(1-a_m);
    nvelo_si(n,:)       = nvelo_si(n,:) + dt * ((1-gamma)*nacc_old(n,:)+gamma*nacc(n,:));
     end
 end
 
% Boundary condition
[nacc]          = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nacc); % Boundary condition for nodal force
[nacc_old]      = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nacc_old); % Boundary condition for nodal force
[nvelo_si]      = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nvelo_si); % Boundary condition for nodal force

%% Update solid particle velocity and position
% [v_ssp,x_sp,d_sp] = Update_Particle_Position(NODES,dt,CONNECT,N,spCount,nmass_si,nforce_si,nmomentum_si,x_spo,v_ssp,x_sp,d_sp);

a_sp(:) = 0;

for sp = 1:spCount
     for j = 1:NODES(sp)
         npid                           = CONNECT{sp}(j);
              if nmass_si(npid)==0
                continue
              end         
         v_ssp(sp,:)    = v_ssp(sp,:) + dt * N{sp}(j) * ((1-gamma)*nacc_old(npid,:)+gamma*nacc(npid,:));
         x_sp(sp,:)     = x_sp(sp,:) + dt * N{sp}(j)*(nvelo_si(npid) + dt * ((0.5-beta)*nacc_old(npid,:)+beta*nacc(npid,:)));
         d_sp(sp,:)     = x_sp(sp,:) - x_spo(sp,:);
         a_sp(sp,:)     = a_sp(sp,:) + N{sp}(j)*nacc(npid,:);
     end   
end
 
% Shape
[~,dN,CONNECT,~,mspoints,NODES] = Compute_Interpolator_CPDI(spCount,cellCount,x_sp,le,NN,LOC,r1_sp,r2_sp,V_sp);

%% Update effective stress
[F_sp,V_sp,s_sp,p_sp] = Update_Stress(CModel,CModel_parameter,...
    NODES,dt,cellCount,mspoints,CONNECT,nvelo_si,dN,...
    F_sp,V_spo,m_sp,s_sp,p_sp,V_sp);

%% Update the topology of particles
[r1_sp,r2_sp] = Update_topology(spCount,F_sp,r1_sp,r10_sp,r2_sp,r20_sp);