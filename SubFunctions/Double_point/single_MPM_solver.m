function[v_ssp,x_sp,F_sp,V_p,s_sp,p_sp,n_sp,k_sp,pore_wp] = single_MPM_solver(CModel,CModel_parameter,nodeCount,spCount,cellCount,...
x_sp,x_spo,d_sp,le,NN,LOC,LOCC,...
b_sp,V_p,ptraction_sp,v_ssp,v_wp,s_sp,pore_wp,m_p,p_sp,k_sp,n_sp,...
nfbcx,nfbcy,fbcx,fbcy,Kw,k,F_sp,V_po,n_o,dt)

%% Store particles into cell
[N,dN,spCONNECT,spElems,mspoints,NODES] = Compute_Interpolator_MPM(spCount,cellCount,x_sp,le,NN,LOC);
% N: Shape function
% dN: Gradient of shape function
% spElems(p): element index where "p" locates
% mspoints(e): all particles indexes where locate in the cell "e"

%% Interpolation from particle to grid task
% Node variables
nmass_si                = zeros(nodeCount,1);                   % Nodal Mass
nmass_wi                = zeros(nodeCount,1);                   % Nodal Mass
nmomentum_si            = zeros(nodeCount,2);                   % Nodal Momentum
nmomentum_wi            = zeros(nodeCount,2);                   % Nodal Momentum
niforce_si              = zeros(nodeCount,2);                   % Nodal Internal force
niforce_wi              = zeros(nodeCount,2);                   % Nodal Internal force
neforce_si              = zeros(nodeCount,2);                   % Nodal External force
neforce_wi              = zeros(nodeCount,2);                   % Nodal External force
traction             = zeros(nodeCount,2);                   % Nodal Traction

% Mass, momentum, external forces
 for p=1:spCount
 % Build stress tensor
 SSP = [s_sp(p,1) s_sp(p,3);s_sp(p,3) s_sp(p,2)];
 
 for j=1:NODES(p)
     npid                           = spCONNECT{p}(j);
     
          if N{p}(j)==1
         continue
          end
 % Mass
 nmass_si(npid)            = nmass_si(npid) + (1-n_sp(p)) * m_sp(p) * N{p}(j);
 nmass_wi(npid)            = nmass_wi(npid) + n_sp(p) * m_wp(p) * N{p}(j);

 % Momentum
 nmomentum_si(npid,:)      = nmomentum_si(npid,:) + (1-n_sp(p)) * m_sp(p) * v_ssp(p,:) * N{p}(j);
 nmomentum_wi(npid,:)      = nmomentum_wi(npid,:) + n_sp(p) * m_wp(p) * v_wp(p,:) * N{p}(j);

 % External forces
 neforce_si(npid,:)         = neforce_si(npid,:) + (1-n_sp(p)) * m_sp(p) * b_sp(p,:) * N{p}(j);
 neforce_wi(npid,:)         = neforce_wi(npid,:) + n_sp(p) * m_wp(p) * b_sp(p,:) * N{p}(j);

 % Traction
 traction(npid,1)       = traction(npid,1) + V_p(p)*ptraction_sp(p,1)*N{p}(j)/le(1);
 traction(npid,2)       = traction(npid,2) + V_p(p)*ptraction_sp(p,2)*N{p}(j)/le(2);
 
  % Internal forces
 niforce_si(npid,:)         = niforce_si(npid,:) - (V_p(p)*SSP*dN{p}(:,j))';
 end 
 end
 
%% Calculate the pore water pressure in cell
active_elements = unique(spElems);
ElemsCount = length (active_elements);

pore_cell            = zeros(ElemsCount,1);
niforce_cell         = zeros(nodeCount,2);

for i=1:ElemsCount
    k = active_elements(i);                 % Index of the element
    Np = length(mspoints{k});               % Number of particles in the element
    for p=1:Np
        pid = mspoints{k}(p);
        pore_cell(i) = pore_cell(i) + pore_wp(pid);
    end
    pore_cell(i) = pore_cell(i)/Np;
end

%% Calculate the internal force of liquid
for i=1:ElemsCount
    k = active_elements(i);                 % Index of the element
    spid    = find(spElems==k);
    % Build pore water pressure tensor
    PORE = [pore_cell(i) 0;0  pore_cell(i)];
    for j=1:NODES(spid(1))
    npid = spCONNECT{spid(1)}(j);
            
            [~,dN1,dN2]=linearshape(LOCC(k,:),LOC(npid,:),le(1),le(2));
          
 % Internal force
 niforce_cell(npid,:) = niforce_cell(npid,:) - le(1)*le(2)*[dN1 dN2]*PORE;
    end
end

%% Calculate the lumped drag matrix Q
drag                 = zeros(nodeCount,1); 
% Calculate the drag force
 for p=1:spCount
 for j=1:NODES(p)
     npid                           = spCONNECT{p}(j);
          if N{p}(j)==0
         continue
          end
% Drag forces
 drag(npid)          = drag(npid) + n_sp(p) * m_p(p)*9.81*N{p}(j)/k_sp(p);
 end 
 end

%% Update momentum
nforce_i = zeros(nodeCount,2);
nacc_si = zeros(nodeCount,2);
nvelo_wi = zeros(nodeCount,2);
% Update total force
nforce_i = niforce_si + niforce_cell +  neforce_i + traction_i;

for n = 1:nodeCount
    if nmass_i(n)==0
        continue
    end
    nacc_si(n,:) = nforce_i(n,:)/nmass_i(n);
end

ndrag                   = zeros(nodeCount,2);
ndrag = niforce_cell + neforce_wi - nmass_wi .* nacc_si;

for n = 1:nodeCount
    if drag(n) ~= 0
        nvelo_wi(n,:) = ndrag(n,:)/drag(n);
    end
end
nmomentum_i   = nmomentum_i + nforce_i*dt;
 
% Boundary condition
[nvelo_wi]  = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nvelo_wi);           % Boundary condition for nodal liquid momentum
[nforce_i]     = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nforce_i);      % Boundary condition for nodal solid force
[nmomentum_i]  = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nmomentum_i);   % Boundary condition for nodal solid momentum

%% Update solid particle velocity and position
[v_ssp,x_sp,~] = Update_Particle_Position(NODES,dt,spCONNECT,N,spCount,nmass_i,nforce_i,nmomentum_i,x_spo,v_ssp,x_sp,d_sp);
% velocity solid particle: v_ssp
% position solid particle: x_sp
  
%% Mapping nodal velocity back to node
% Solid phase
[nvelo_i] = Interpolate_velocity_back(NODES,nodeCount,spCount,spCONNECT,m_p,v_ssp,N,nmass_i);
% Boundary condition
[nvelo_i] = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nvelo_i); % Boundary condition for nodal force

%% Update effective stress
[F_sp,V_p,s_sp,p_sp] = Update_Stress(CModel,CModel_parameter,...
    NODES,dt,cellCount,mspoints,spCONNECT,nvelo_i,dN,...
    F_sp,V_po,m_p,s_sp,p_sp,V_p);

%% Update pore water pressure
[pore_wp] = Update_Pore_Pressure(spCount,NODES,dt,cellCount,mspoints,spCONNECT,nvelo_i,nvelo_wi,n_sp,pore_wp,dN,Kw);

%% Update porousity
[n_sp,k_sp]=Update_porosity(n_o,F_sp,k,spCount);
