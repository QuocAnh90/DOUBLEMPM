function[N_MLS,singular_node]=CubicNode_IMLS(spCount,nodeCount,npoints,CONNECT_TEMP_N,x_sp,LOC,le)

% Funtion for mapping from particles to nodes using cubic spline
Nn_local        = zeros(spCount,16);              % Value of shape function

singular_node = [];

for n = 1:nodeCount
     if isempty(npoints{n})==1
         continue
     end
     cv = x_sp(npoints{n},:);
     cv = cv';
     Np = length(npoints{n});
     won = ones(1,Np);
     wn = zeros(1,Np);
     for i = 1:Np
         pid = npoints{n}(i);
%          N_local(pid,n) = linearshape(x_sp(pid,1:2),LOC(n,:),le(1),le(2));
         Nn_local(pid,n) = Cubic_Bspline(x_sp(pid,1:2),LOC(n,:),le(1),le(2));
         wn(i) = Nn_local(pid,n);
     end
     p = [won;cv];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q=zeros(3,Np);
q(1,1:Np)=1;
gama=zeros(1,3);
bita=zeros(1,3);
Rfa=zeros(1,3);

p1xi=p(1,1:Np);
q1xi=q(1,1:Np);

p2xi=p(2,1:Np);
gama1=0;
for i=1:Np
    gama1=gama1+wn(i)*q1xi(i)^2;
end
gama(1)=gama1;
bita21=0;
for i=1:Np
    bita21=bita21+wn(i)*p2xi(i)*q1xi(i);
end
bita(1)=bita21;
if gama1==0
    Rfa21=0;
    singular_node = [singular_node n];
else
    Rfa21=bita21/gama1;
end
Rfa(1)=Rfa21;

for i=1:Np
    q(2,i)=p2xi(i)-Rfa21*q1xi(i);
end
q2xi=q(2,1:Np);

p3xi=p(3,1:Np);
gama2=0;
bita31=0;
bita32=0;
for i=1:Np
    gama2=gama2+wn(i)*q2xi(i)^2;
    bita31=bita31+wn(i)*p3xi(i)*q1xi(i);
    bita32=bita32+wn(i)*p3xi(i)*q2xi(i);
end
gama(2)=gama2;
bita(2)=bita31;
bita(3)=bita32;
if gama1==0
    Rfa31=0;
    singular_node = [singular_node n];
else
    Rfa31=bita31/gama1;
end

if gama2==0
    Rfa32=0;
    singular_node = [singular_node n];
else
    Rfa32=bita32/gama2;
end

Rfa(2)=Rfa31;
Rfa(3)=Rfa32;
for i=1:Np
    q(3,i)=p3xi(i)-Rfa31*q1xi(i)-Rfa32*q2xi(i);
end

q3xi=q(3,1:Np);

gama3=0;
for i=1:Np
    gama3=gama3+wn(i)*q3xi(i)^2;
end
gama(3)=gama3;

gposx=LOC(n,1);
gposy=LOC(n,2);
qq=zeros(1,3);
qq(1)=1;
qq(2)=gposx-Rfa(1);
qq(3)=gposy-(Rfa(2)*qq(1)+Rfa(3)*qq(2));
cc=zeros(3,Np);
for j=1:3
    for i=1:Np
        if gama(j)==0
        cc(j,i)=0;
        else
        cc(j,i)=(qq(j)*q(j,i))/gama(j);
        end
    end
end
phi=zeros(1,Np);
for i=1:Np
    cji=0;
    for j=1:3
        cji=cji+cc(j,i);
    end
    phi(i)=wn(i)*cji;
end

     for i = 1:Np
         pid = npoints{n}(i);
%      id_n = find(CONNECT_TEMP(pid,:)==n);
     id_n = find(CONNECT_TEMP_N(pid,:)==n);
     N_MLS{pid}(id_n) = phi(i);
     end
 end