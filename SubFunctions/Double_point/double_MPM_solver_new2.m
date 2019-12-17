function[v_ssp,x_sp,F_sp,V_sp,s_sp,p_sp,v_wp,x_wp,pore_wp,F_wp,V_wp] = double_MPM_solver_new2(CModel,CModel_parameter,nodeCount,wpCount,spCount,cellCount,...
x_wp,x_sp,x_wpo,x_spo,d_wp,d_sp,F_wp,le,NN,LOC,LOCC,...
b_sp,b_wp,V_sp,ptraction_sp,ptraction_wp,v_ssp,v_wp,s_sp,pore_wp,m_wp,m_sp,p_sp,...
nfbcx,nfbcy,fbcx,fbcy,k,n_o,F_sp,V_spo,V_wp,dt,psp)

%% Store particles into cell
% Solid
[Ns,dNs,spCONNECT,spElems,mspoints,NODES_sp] = Compute_Interpolator_MPM(spCount,cellCount,x_sp,le,NN,LOC);

active_elements_solid = unique(spElems);
ElemsCount_solid = length (active_elements_solid);
% Liquid
[Nw,dNw,wpCONNECT,wpElems,mwpoints,NODES_wp] = Compute_Interpolator_MPM(wpCount,cellCount,x_wp,le,NN,LOC);
active_elements_liquid = unique(wpElems);
ElemsCount_liquid = length (active_elements_liquid);
% N: Shape function
% dN: Gradient of shape function
% spElems(p): element index where "p" locates
% mspoints(e): all particles indexes where locate in the cell "e"

%% Mapping from particle to nodes
% Solid variables
[nmass_si,nmomentum_si,niforce_si,neforce_si,traction_si,nvolume_si]=Interpolate_Particle_To_Grid(...
 NODES_sp,nodeCount,spCONNECT,le,Ns,dNs,spCount,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp);
 
% Liquid variables
[nmass_wi,nmomentum_wi,niforce_wi,neforce_wi,traction_wi]=Interpolate_Liquid_Particle_To_Grid_new2(...
NODES_wp,nodeCount,wpCONNECT,le,Nw,dNw,wpCount,wpElems,mwpoints,LOCC,LOC,b_wp,V_wp,ptraction_wp,v_wp,pore_wp,m_wp,active_elements_liquid,ElemsCount_liquid);

%% Update nodal Porosity and permeability
% Node variables
n_i                 = ones(nodeCount,1);                   % Nodal porousity
k_i                 = zeros(nodeCount,1);                   % Nodal permeability
ndrag               = zeros(nodeCount,2);
nvelo_si            = zeros(nodeCount,2);
nvelo_wi            = zeros(nodeCount,2);

 for n = 1:nodeCount 
 % Porosity
 if nvolume_si(n) ~=0
 n_i(n)            = 1-(nmass_si(n))/(psp*nvolume_si(n));
 end
 
 % Momentum
 k_i(n)            = k;
 % Solid velocity
 if nmass_si(n) ~=0
 nvelo_si(n,:)     = nmomentum_si(n,:)/nmass_si(n);
 end
 % Liquid velocity
 if nmass_wi(n) ~=0
 nvelo_wi(n,:)     = nmomentum_wi(n,:)/nmass_wi(n);
 end
 % Drag forces
 if nvolume_si(n) ~=0
 ndrag(n,:)        = nmass_wi(n)*9.81/k_i(n)*(nvelo_wi(n,:) - nvelo_si(n,:));
 end 
 end

 %% Update momentum
 nacc_si              = zeros(nodeCount,2);
 nacc_wi              = zeros(nodeCount,2);
 
 % Update force and momentum for liquid
 nforce_wi =  n_i .* niforce_wi + neforce_wi + traction_wi - n_i.*ndrag;

 for n = 1:nodeCount
 if nmass_wi(n) ~=0
 nacc_wi(n,:) =  nforce_wi(n,:)/nmass_wi(n);
 end
 end
 nvelo_wi = nvelo_wi + nacc_wi*dt;

 % Update force and momentum for liqid
 for n = 1:nodeCount
%  if nvolume_si(n) ~=0
 nforce_si(n,:) = niforce_si(n,:) + niforce_wi(n,:) +  neforce_si(n,:) + neforce_wi(n,:) + traction_wi(n,:) + traction_si(n,:) - nforce_wi(n,:);
%  end
 end

 for n = 1:nodeCount
 if nmass_si(n) ~=0
 nacc_si(n,:) =  nforce_si(n,:)/nmass_si(n);
 end
 end
 nvelo_si   = nvelo_si + nacc_si*dt;
 
 % Boundary condition
 [nacc_wi]     = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nacc_wi);      % Boundary condition for nodal liquid force
 [nacc_si]     = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nacc_si);      % Boundary condition for nodal solid force
 [nvelo_si]      = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nvelo_si);
 [nvelo_wi]      = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nvelo_wi);

%% Update solid particle velocity and position
[v_ssp,x_sp] = Update_Particle_Position_new(NODES_sp,dt,spCONNECT,Ns,spCount,nmass_si,nacc_si,nvelo_si,n_i,x_spo,v_ssp,x_sp,d_sp);
% velocity solid particle: v_ssp
% position solid particle: x_sp

[v_wp,x_wp,n_wp] = Update_Particle_Position_new(NODES_wp,dt,wpCONNECT,Nw,wpCount,nmass_wi,nacc_wi,nvelo_wi,n_i,x_wpo,v_wp,x_wp,d_wp);
% velocity liquid particle: v_ssp
% position liquid particle: x_sp

%% Update effective stress
[F_sp,V_sp,s_sp,p_sp] = Update_Stress(CModel,CModel_parameter,...
    NODES_sp,dt,active_elements_solid,ElemsCount_solid,mspoints,spCONNECT,nvelo_si,dNs,...
    F_sp,V_spo,m_sp,s_sp,p_sp,V_sp);

%% Update pore water pressure
% Kw                      = 2200e6         ;
% [pore_wp] = Update_Pore_Pressure_new1(...
%     wpCount,NODES_wp,dt,active_elements_liquid,ElemsCount_liquid,nodeCount,mwpoints,wpCONNECT,nvelo_si,nvelo_wi,n_wp,n_i,pore_wp,dNw,Kw);

[F_wp,pore_wp] = Update_Pore_Pressure_new2(...
    wpCount,NODES_wp,dt,active_elements_liquid,ElemsCount_liquid,nodeCount,mwpoints,wpCONNECT,nvelo_si,nvelo_wi,n_i,n_wp,n_o,pore_wp,F_wp,dNw);
test = 1;