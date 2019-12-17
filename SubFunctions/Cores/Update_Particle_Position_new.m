function [v_ssp,x_sp,n_sp,d_sp] = Update_Particle_Position_new(NODES,dt,CONNECT,N,spCount,nmass,nacc,nvelo,n_i,x_spo,v_ssp,x_sp,d_sp)

n_sp                    = zeros(spCount,1);                 % liquid particle porosity

for p = 1:spCount
     for j = 1:NODES(p)
         npid                           = CONNECT{p}(j);
              if nmass(npid)==0
                continue
              end         
         v_ssp(p,:)                      = v_ssp(p,:) + dt * N{p}(j) * nacc(npid,:);
         x_sp(p,:)                       = x_sp(p,:) + nvelo(npid,:)*N{p}(j)*dt;
         d_sp(p,:)                       = x_sp(p,:) - x_spo(p,:);
         n_sp(p)                         = n_sp(p) + n_i(npid) * N{p}(j);
     end   
end