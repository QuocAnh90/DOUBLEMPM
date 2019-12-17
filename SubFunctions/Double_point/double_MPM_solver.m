function[v_ssp,x_sp,F_sp,V_sp,s_sp,p_sp,n_sp,k_sp,v_wp,x_wp,pore_wp] = double_MPM_solver(CModel,CModel_parameter,nodeCount,wpCount,spCount,cellCount,...
x_wp,x_sp,x_wpo,x_spo,d_wp,d_sp,le,NN,LOC,LOCC,...
b_sp,b_wp,V_sp,ptraction_sp,ptraction_wp,v_ssp,v_wp,s_sp,pore_wp,m_wp,m_sp,p_sp,k_sp,n_sp,...
nfbcx,nfbcy,fbcx,fbcy,Kw,k,F_sp,V_spo,V_wpo,n_o,dt)

%% Store particles into cell
% Solid
[Ns,dNs,spCONNECT,spElems,mspoints,NODES_sp] = Compute_Interpolator_MPM(spCount,cellCount,x_sp,le,NN,LOC);
% Liquid
[Nw,dNw,wpCONNECT,wpElems,mwpoints,NODES_wp] = Compute_Interpolator_MPM(wpCount,cellCount,x_wp,le,NN,LOC);
% N: Shape function
% dN: Gradient of shape function
% spElems(p): element index where "p" locates
% mspoints(e): all particles indexes where locate in the cell "e"
active_elements_solid = unique(spElems);
ElemsCount_solid = length (active_elements_solid);
active_elements_liquid = unique(wpElems);
ElemsCount_liquid = length (active_elements_liquid);

 %% Mapping from particle to nodes
 % Solid variables
[nmass_si,nmomentum_si,niforce_si,neforce_si,traction_si,nvolume_si]=Interpolate_Particle_To_Grid(NODES_sp,nodeCount,spCONNECT,le,Ns,dNs,spCount,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp);

% Update nodal Porousity and permeability
[n_si,dn_si,k_si]=Interpolate_porosity_to_grid(NODES_sp,nodeCount,spCONNECT,dNs,Ns,spCount,k_sp,n_sp,m_sp,nmass_si);

% Update liquid particle porosity, permeability and volumn
[n_wp,k_wp,V_wp]=Update_liquid_porosity_particle(NODES_wp,wpCONNECT,Nw,wpCount,k_si,n_si,nmass_si,n_o,V_wpo);

% Liquid variables
[nmass_wi,nmomentum_wi,niforce_cell,neforce_wi,traction_wi,alpha]=Interpolate_Liquid_Particle_To_Grid(NODES_wp,nodeCount,wpCONNECT,mwpoints,wpElems,LOCC,LOC,le,Nw,wpCount,b_wp,V_wp,ptraction_wp,v_wp,pore_wp,m_wp,n_wp);

% Calculate the internal force
[ndrag,niforce_wi] = Calculate_drag_force(NODES_wp,nodeCount,wpCONNECT,Nw,wpCount,m_wp,k_wp,nmass_wi,nmass_si,nmomentum_wi,nmomentum_si,alpha,niforce_cell);

%% Update momentum
% Update force and momentum for liqid
for n=1:nodeCount
nforce_wi(n,:) =  niforce_wi(n,:) + neforce_wi(n,:) + traction_wi(n,:) - alpha(n)*ndrag(n,:);
end

nmomentum_wi = nmomentum_wi + nforce_wi*dt;

% nforce_si = niforce_si + niforce_cell +  neforce_si + neforce_wi + traction_wi + traction_si - nforce_wi;
nforce_si = niforce_si + (1-alpha).*niforce_cell +  neforce_si + traction_si + alpha.*ndrag;

nmomentum_si   = nmomentum_si + nforce_si*dt;
 
nvelo_wi                = zeros(nodeCount,2); 
nvelo_si                = zeros(nodeCount,2); 

for n = 1:nodeCount
    if nmass_wi(n)==0
        continue
    end
    nvelo_wi(n,:) = nmomentum_wi(n,:)/nmass_wi(n);
end

for n = 1:nodeCount
    if nmass_si(n)==0
        continue
    end
nvelo_si(n,:) = nmomentum_si(n,:)/nmass_si(n);
end

% Boundary condition
[nforce_wi]     = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nforce_wi);      % Boundary condition for nodal liquid force
[nmomentum_wi]  = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nmomentum_wi);   % Boundary condition for nodal liquid momentum
[nforce_si]     = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nforce_si);      % Boundary condition for nodal solid force
[nmomentum_si]  = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nmomentum_si);   % Boundary condition for nodal solid momentum
[nvelo_si]      = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nvelo_si);
[nvelo_wi]      = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nvelo_wi);

%% Update solid particle velocity and position
[v_ssp,x_sp,~] = Update_Particle_Position(NODES_sp,dt,spCONNECT,Ns,spCount,nmass_si,nforce_si,nmomentum_si,x_spo,v_ssp,x_sp,d_sp);
% velocity solid particle: v_ssp
% position solid particle: x_sp

[v_wp,x_wp,~] = Update_Particle_Position(NODES_wp,dt,wpCONNECT,Nw,wpCount,nmass_wi,nforce_wi,nmomentum_wi,x_wpo,v_wp,x_wp,d_wp);
% velocity liquid particle: v_ssp
% position liquid particle: x_sp
  
% %% Mapping nodal velocity back to node
% % Solid phase
% [nvelo_si] = Interpolate_velocity_back(NODES_sp,nodeCount,spCount,spCONNECT,m_sp,v_ssp,Ns,nmass_si);
% % Boundary condition
% [nvelo_si] = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nvelo_si); % Boundary condition for nodal force
% 
% % Liquid phase
% [nvelo_wi] = Interpolate_velocity_back(NODES_wp,nodeCount,wpCount,wpCONNECT,m_wp,v_wp,Nw,nmass_wi);
% 
% % Boundary condition
% [nvelo_wi] = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nvelo_wi); % Boundary condition for nodal force

%% Update effective stress
[F_sp,V_sp,s_sp,p_sp] = Update_Stress(CModel,CModel_parameter,...
    NODES_sp,dt,active_elements_solid,ElemsCount_solid,mspoints,spCONNECT,nvelo_si,dNs,...
    F_sp,V_spo,m_sp,s_sp,p_sp,V_sp);

[n_sp,k_sp]=Update_porosity(n_o,F_sp,k,spCount);

%% Update pore water pressure
[pore_wp] = Update_Pore_Pressure(wpCount,NODES_wp,dt,cellCount,mwpoints,wpCONNECT,nvelo_si,nvelo_wi,n_wp,pore_wp,dNw,Kw);

% [pore_wp] = Update_Pore_Pressure_Filter(wpCount,NODES_wp,dt,cellCount,mwpoints,wpCONNECT,nvelo_si,nvelo_wi,n_wp,pore_wp,dNw,Kw);
