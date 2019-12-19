function [nmass_wi,nmomentum_wi,niforce_wi,neforce_wi,traction_wi]=Interpolate_Liquid_Particle_To_Grid(...
NODES_wp,nodeCount,wpCONNECT,le,Nw,dNw,wpCount,b_wp,V_wp,ptraction_wp,v_wp,pore_wp,m_wp)
%% Interpolation from particle to grid task
% Node variables
nmass_wi                = zeros(nodeCount,1);                   % Nodal Mass
nmomentum_wi            = zeros(nodeCount,2);                   % Nodal Momentum
neforce_wi              = zeros(nodeCount,2);                   % Nodal External force
traction_wi             = zeros(nodeCount,2);                   % Nodal Traction
niforce_wi              = zeros(nodeCount,2);

 for p=1:wpCount
     
 % Build stress tensor
%  PORE = [pore_wp(p) 0;0 pore_wp(p)];
 PORE = [pore_wp(p,1) pore_wp(p,3);pore_wp(p,4) pore_wp(p,2)];
 
 for j=1:NODES_wp(p)
     npid                           = wpCONNECT{p}(j);
     
          if Nw{p}(j)==0
         continue
          end
     
 % Mass
 nmass_wi(npid)            = nmass_wi(npid) + m_wp(p)*Nw{p}(j);
 
 % Momentumlp
 nmomentum_wi(npid,:)      = nmomentum_wi(npid,:) + m_wp(p)*v_wp(p,:)*Nw{p}(j);

 % Internal forces
 niforce_wi(npid,:)         = niforce_wi(npid,:) - (V_wp(p)*PORE*dNw{p}(:,j))';
 
 % External forces
 neforce_wi(npid,:)         = neforce_wi(npid,:) + b_wp(p,:)*m_wp(p)*Nw{p}(j);

 % Traction
 traction_wi(npid,1)       = traction_wi(npid,1) + V_wp(p)*ptraction_wp(p,1)*Nw{p}(j)/le(2);
 traction_wi(npid,2)       = traction_wi(npid,2) + V_wp(p)*ptraction_wp(p,2)*Nw{p}(j)/le(1); 
 end 
 end