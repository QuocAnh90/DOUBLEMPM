function [nmass_si,nmomentum_si,niforce_si,neforce_si,traction_si,nvolume_si]=Interpolate_Particle_To_Grid(...
    NODES_sp,nodeCount,spCONNECT,le,Ns,dNs,spCount,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp)

%% Mapping from solid particle to nodes
% Node variables
nmass_si                = zeros(nodeCount,1);                   % Nodal Mass
nvolume_si                = zeros(nodeCount,1);                   % Nodal Mass
nmomentum_si            = zeros(nodeCount,2);                   % Nodal Momentum
niforce_si              = zeros(nodeCount,2);                   % Nodal Internal force
neforce_si              = zeros(nodeCount,2);                   % Nodal External force
traction_si             = zeros(nodeCount,2);                   % Nodal Traction

 for p=1:spCount
 % Build stress tensor
 SSP = [s_sp(p,1) s_sp(p,3);s_sp(p,3) s_sp(p,2)];
 
 for j=1:NODES_sp(p)
     npid                           = spCONNECT{p}(j);
     
          if Ns{p}(j)==1
         continue
          end

 % Mass
 nmass_si(npid)            = nmass_si(npid) + m_sp(p)*Ns{p}(j);
 
 % Mass
 nvolume_si(npid)            = nvolume_si(npid) + V_sp(p)*Ns{p}(j);
 
 % Momentum
 nmomentum_si(npid,:)      = nmomentum_si(npid,:) + m_sp(p)*v_ssp(p,:)*Ns{p}(j);
 
 % Internal forces
 niforce_si(npid,:)         = niforce_si(npid,:) - (V_sp(p)*SSP*dNs{p}(:,j))';

 % External forces
 neforce_si(npid,:)         = neforce_si(npid,:) + b_sp(p,:)*m_sp(p)*Ns{p}(j);

 % Traction
 traction_si(npid,1)       = traction_si(npid,1) - V_sp(p)*ptraction_sp(p,1)*Ns{p}(j)/le(1);
 traction_si(npid,2)       = traction_si(npid,2) - V_sp(p)*ptraction_sp(p,2)*Ns{p}(j)/le(2); 
 end 
 end
 
 % % Notes on traction calculation
% % If traction / le: it means that the traction layer is equal to the cell
% % layer
% % If traction / lp: it means that the traction layer is equal to the
% % particle  layer