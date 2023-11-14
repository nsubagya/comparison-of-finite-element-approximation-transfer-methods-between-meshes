%x-FE space 1
%y-FE space 2
%u-sulution set in x
%graphflag- if graphflag >0 graph is on O/W off
function u2=L2proj7(x,y,u)
%Example command to type in command window
%L2proj5(0:0.01:10,0:0.01:10,1) should give us graphs top of each other
z=union(x,y);
%end the process if first and last nodes are not equl in x, y, z
u1=zeros(1,length(z));
%Initiate mass matrix and initiate the projection matrix for FE space 2
M=sparse(1:length(y),1:length(y),0,length(y),length(y));
b=zeros(1,length(y));
for i=1:length(y)-1
    H(i)=y(i+1)-y(i);
    M(i:i+1,i:i+1)=M(i:i+1,i:i+1)+(H(i)/6)*[2 1;1 2];
end

%initiate the local intervals I=[i j k] of k interval
I=[0 0 0];
for k=1:length(z)-1
    switch k
        case 1
            i=1;j=1;
        otherwise
            if z(k)>=x(i+1)
                i=i+1;
            end
            if z(k)>=y(j+1)
                j=j+1;
            end
    end
    %I=[i j k] (keeping this to troubleshooting if needed)
    %ratio for interpolate (z(k) and z(k+1)) LOCALLY
    lambdax=(z(k)-x(i))/(x(i+1)-x(i));
    alphax=(z(k+1)-x(i))/(x(i+1)-x(i));
    u1(k)=lambdax*u(i+1)+(1-lambdax)*u(i);
    u1(k+1)=alphax*u(i+1)+(1-alphax)*u(i);
    
    %ratio for interpolate two (phi(k) and phi(k+1)) LOCALLY
    shyy=(z(k)-y(j))/(y(j+1)-y(j));
    phiy=(z(k+1)-y(j))/(y(j+1)-y(j));
    
    %both hat funcion values for local k node
    phi1(1)=(1-shyy);
    phi1(2)=(1-phiy);
    
    phi2(1)=shyy;
    phi2(2)=phiy;
    
    %find RHS values
    b(j)=b(j)+(u1(k)*phi1(1)+4*(u1(k)+u1(k+1))/2 *(phi1(1)+phi1(2))/2+u1(k+1)*phi1(2))*(z(k+1)-z(k))/6;
    b(j+1)=b(j+1)+(u1(k)*phi2(1)+4*(u1(k)+u1(k+1))/2 *(phi2(1)+phi2(2))/2+u1(k+1)*phi2(2))*(z(k+1)-z(k))/6;
end

u2=b/M;
end