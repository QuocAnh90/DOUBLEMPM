function[v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp,r1_sp,r2_sp] = CPDI_solver_PIC(CModel,CModel_parameter,nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt,r1_sp,r10_sp,r2_sp,r20_sp,t)
     
%% Store particles into cell
[N,dN,CONNECT,spElems,mspoints,NODES] = Compute_Interpolator_CPDI(spCount,cellCount,x_sp,le,NN,LOC,r1_sp,r2_sp,V_sp);
% N: Shape function
% dN: Gradient of shape function
% spElems(p): element index where "p" locates
% mspoints(e): all particles indexes where locate in the cell "e"

 %% Mapping from particle to nodes
[nmass_si,nmomentum_si,niforce_si,neforce_si,traction_si]=Interpolate_Particle_To_Grid(NODES,nodeCount,CONNECT,le,N,dN,spCount,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp);
% nmass_si: nodal mass
% nmomentum_si: nodal momentum
% niforce_si: nodal internal force
% neforce_si: nodal external force
% traction_si: nodal traction

%% Update momentum
% Update force and momentum
 nforce_si      = niforce_si + neforce_si + traction_si;
 nmomentum_si   = nmomentum_si + nforce_si*dt;
 
nvelo_si = zeros(nodeCount,2);
for n = 1:nodeCount
    if nmass_si(n) >1e-12
    nvelo_si(n,:) = nmomentum_si(n,:)/nmass_si(n);
    end
end
 
 
% Boundary condition
[nforce_si]     = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nforce_si); % Boundary condition for nodal force
% [nmomentum_si]  = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nmomentum_si); % Boundary condition for nodal force
[nvelo_si] = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nvelo_si); % Boundary condition for nodal force

% for n=1:nodeCount
%     if LOC(n,1)<2.5*le(1)
%         if t<0.0000012241
%         nvelo_si(n,1) = 0.00002;
%         end
%     end
% end

%% Update solid particle velocity and position
[~,x_sp,d_sp] = Update_Particle_Positionnull(NODES,dt,CONNECT,N,spCount,nmass_si,nforce_si,nvelo_si,x_spo,v_ssp,x_sp,d_sp);
% velocity particle: v_ssp
% position particle: x_sp
% displacement particle: d_sp
 

% %% PIC filter
% %% Interpolation from particle to grid task
% % Node variables
% nmomentum            = zeros(nodeCount,2);                   % Nodal Momentum
%  for p=1:spCount 
%  for j=1:NODES(p)
%      npid                           = CONNECT{p}(j);
%      
%           if N{p}(j)==1
%          continue
%           end
%      
%  % Momentum
%  nmomentum(npid,:)      = nmomentum(npid,:) + m_sp(p)*v_ssp(p,:)*N{p}(j);
%  end 
%  end
%  
%  nvelo_si11 = zeros(nodeCount,2);
% for n = 1:nodeCount
%     if nmass_si(n) >1e-12
%     nvelo_si11(n,:) = nmomentum_si(n,:)/nmass_si(n);
%     end
% end

v_ssp = zeros(spCount,2);

for sp = 1:spCount
     for j = 1:NODES(sp)
         npid                           = CONNECT{sp}(j);
              if nmass_si(npid)==0
                continue
              end         
         v_ssp(sp,:)                      = v_ssp(sp,:) + N{sp}(j) * nvelo_si(npid,:);
     end   
end

 
%% Mapping nodal velocity back to node
% [nvelo_si] = Interpolate_velocity_back(NODES,nodeCount,spCount,CONNECT,m_sp,v_ssp,N,nmass_si);
% nvelo_si = zeros(nodeCount,2);
% for n = 1:nodeCount
%     if nmass_si(n) >1e-12
%     nvelo_si(n,:) = nmomentum_si(n,:)/nmass_si(n);
%     end
% end

% Boundary condition
[~,dN,CONNECT,~,mspoints,NODES] = Compute_Interpolator_CPDI(spCount,cellCount,x_sp,le,NN,LOC,r1_sp,r2_sp,V_sp);

%% Update effective stress
[F_sp,V_sp,s_sp,p_sp] = Update_Stress(CModel,CModel_parameter,...
    NODES,dt,cellCount,mspoints,CONNECT,nvelo_si,dN,...
    F_sp,V_spo,m_sp,s_sp,p_sp,V_sp);

%% Update the topology of particles
[r1_sp,r2_sp] = Update_topology(spCount,F_sp,r1_sp,r10_sp,r2_sp,r20_sp);

