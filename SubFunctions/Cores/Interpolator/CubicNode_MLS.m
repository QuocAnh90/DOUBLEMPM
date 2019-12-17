function N_MLS = CubicNode_MLS(nodeCount,spCount,NODES_N,CONNECT_TEMP_N,CONNECT_N,x_sp,LOC,le,r,r_n)

% Funtion for mapping from particles to nodes using linear basis
Nn_local         = zeros(spCount,16);              % Value of shape function

% Funtion for mapping from particles to nodes using moving least square
% Least square matrix
A               = cell(nodeCount,1);

for i = 1:nodeCount
    A{i} = zeros(r_n(i),r_n(i));
end

for p = 1:spCount
for i = 1:NODES_N(p)
    npid = CONNECT_TEMP_N(p,i);
     % Compute the shape functions and gradient of the shape functions
    [Nn_local(p,i)]=Cubic_Bspline(x_sp(p,1:2),LOC(npid,:),le(1),le(2));
end
 end

% Compute A
for p = 1:spCount
     if r==1
         Pxy = 1;
     elseif r==3
         Pxy = [1 ; x_sp(p,1) ; x_sp(p,2)];    
     elseif r==6
         Pxy = [1 ; x_sp(p,1) ; x_sp(p,2) ; x_sp(p,1)^2 ; x_sp(p,1)*x_sp(p,2) ; x_sp(p,2)^2];   
     end
         
     PPxy = Pxy * Pxy';
     
     for i = 1:NODES_N(p)
         npid = CONNECT_N{p}(i);
         A{npid} = A{npid} + Nn_local(p,i) * PPxy;  
     end
end

% % Check rank D
% for n=1:nodeCount
% %     if rank(A{n})<6
%     if rcond(A{n})<=9.999999e-10
%         r_n(n) = 3;
%     end
% end

% Compute shape function and gradients
 for p = 1:spCount
     if r==1
         Pxy = 1;
     elseif r==3
         Pxy = [1 ; x_sp(p,1) ; x_sp(p,2)];   
     elseif r==6
         Pxy = [1 ; x_sp(p,1) ; x_sp(p,2) ; x_sp(p,1)^2 ; x_sp(p,1)*x_sp(p,2) ; x_sp(p,2)^2];   
     end
     
     for i = 1:NODES_N(p)
         npid = CONNECT_N{p}(i);
         if r==1
         pn = 1;
         elseif r==3
         pn = [1 ; LOC(npid,1) ; LOC(npid,2)];   
         elseif r==6
         pn = [1 ; LOC(npid,1) ; LOC(npid,2) ; LOC(npid,1)^2 ; LOC(npid,1)*LOC(npid,1) ; LOC(npid,1)^2];   
         end
         
         if r_n(npid)==6
         rn = A{npid} \ pn;
         N_MLS{p}(i) = rn' * Nn_local(p,i) * Pxy;
         
         elseif r_n(npid)==3
         rn = A{npid}(1:3,1:3) \ pn(1:3);
         N_MLS{p}(i) = rn' * Nn_local(p,i) * Pxy(1:3);
           
         end
         
     end
 end