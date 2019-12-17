function [n_si,dn_si,k_si]=Interpolate_porosity_to_grid(NODES,nodeCount,CONNECT,dN,N,spCount,k_sp,n_sp,m_sp,nmass_si)

%% Interpolation from particle to grid task
% Node variables
n_si                = zeros(nodeCount,1);                   % Nodal porousity
k_si                = zeros(nodeCount,1);                   % Nodal permeability
dn_si                = zeros(nodeCount,2); 

 for p=1:spCount
 for j=1:NODES(p) 
     npid                           = CONNECT{p}(j);
     
          if N{p}(j)==0
         continue
          end
     
 % Porosity
 n_si(npid)            = n_si(npid) +  n_sp(p)*m_sp(p)*N{p}(j);
 
 % Porosity gradient
%  dn_si(npid,:)            = dn_si(npid,:) +  n_sp(p)*m_sp(p)*dN{p}(:,j)';
 
 % Momentum
 k_si(npid)            = k_si(npid) +  k_sp(p)*m_sp(p)*N{p}(j);
 end 
 end
 
 for n=1:nodeCount   
     if nmass_si(n)==0
         continue
     end
    n_si(n) = n_si(n)/nmass_si(n);
%     dn_si(n,:) = dn_si(n,:)/nmass_si(n);
    k_si(n) = k_si(n)/nmass_si(n);
end