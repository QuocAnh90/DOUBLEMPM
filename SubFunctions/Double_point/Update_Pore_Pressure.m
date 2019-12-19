function [F_wp,pore_wp] = Update_Pore_Pressure_new2(...
    wpCount,NODES,dt,active_elements,ElemsCount,nodeCount,mwpoints,CONNECT,nvelo_si,nvelo_wi,n_i,n_wp,n_o,pore_wp,F_wp,dN)

de_wp = zeros(nodeCount,1);

for n = 1:nodeCount
 if norm(nvelo_si(n,:)) ~=0
nvelo_mixture(n,:) = (1-n_i(n))*nvelo_wi(n,:) - n_i(n)*nvelo_si(n,:); 
 else
 nvelo_mixture = nvelo_wi;
 end
end
 
L_mixture = cell(wpCount,1);
dWP   = cell(wpCount,1); 

% Calculate pore water pressure for liquid phase
for i=1:ElemsCount
    c = active_elements(i);                                                 % loop all the element
    mptw = mwpoints{c};                                                % all particles in the element
    
    for wp = 1:length(mptw)                                            % loop all particle in that element
        wpid = mptw(wp);                                               % wpid: global index of the particle
        L_mixture{wpid} = zeros(2,2);                                  % Strain rate tensor of liquid phase
    
for j=1:NODES(wpid)                                                    % loop all nodes interacting with particle wpid
          if dN{wpid}(j)==0
         continue
          end
            npid = CONNECT{wpid}(j);                                   % npid: global index of node
            L_mixture{wpid} = L_mixture{wpid} + (dN{wpid}(:,j)*nvelo_mixture(npid,:));          % Interpolate the strain rate from nodes to particles
end
        
        dWP{wpid} = dt*(L_mixture{wpid} + L_mixture{wpid}')/2;      
        de_wp(wpid)  = 1/n_wp(wpid) * trace(dWP{wpid});               
%         pore_wp(wpid,1) = pore_wp(wpid,1) + 220000000 * de_wp(wpid);
%         pore_wp(wpid,2) = pore_wp(wpid,2) + 220000000 * de_wp(wpid);
        
        % Input
        K                       = 150000;                          % Bulk modulus
        u                       = 0.5;                            % vicosity
        gamma                   = 7.0;                            % gamma
        
         
        % Rate of deformation tensor
        D = (L_mixture{wpid} + L_mixture{wpid}')/2;
%         J = (1-n_o)/(1-n_wp(wpid));
        F_wp{wpid} = (eye(2,2)+L_mixture{wpid}*dt)*F_wp{wpid};                           
        J = det(F_wp{wpid});
        
        % The Deviatoric part of the rate of deformation tensor
        Dprime = D - 1/3*trace(D)*eye(2,2);
        Shear = 2*u*Dprime;

        % The isotropic pressure
        p = K*(J^gamma-1);

        % Cauchy stress 
        PORE = p*eye(2,2) + Shear;

        % Regroup the stress
        pore_wp(wpid,1) = PORE(1,1);
        pore_wp(wpid,2) = PORE(2,2);
        pore_wp(wpid,3) = PORE(1,2);
        pore_wp(wpid,4) = PORE(2,1);
    end      
end