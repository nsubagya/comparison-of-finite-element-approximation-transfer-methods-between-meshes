%x-mesh x
%y-mesh y
%u-initial function on mesh x
%we need Lagrangeint.m and gauss.m
function u2=L2projwosupmesh8(x,y,u,gausspt)
%initiate mass matrix and RHS.
b=zeros(length(y),1);
M=sparse(1:length(y),1:length(y),0,length(y),length(y));

%mass matrix for mesh y
for i=1:length(y)-1
    H=y(i+1)-y(i);
    M(i:i+1,i:i+1)=M(i:i+1,i:i+1)+(H/6)*[2 1;1 2];
end
%building the RHS on mesh y
for j=1:length(y)-1
    %import integration (u*phi) of each element of y from gauss.m
   u3=gauss(x,y,j,u,gausspt);
   b(j)=b(j)+u3(1);
   b(j+1)=b(j+1)+u3(2);
end
%solve M*u2=b system by using \
u2=M\b;
u2=u2';%just to use it in heat and transport equation