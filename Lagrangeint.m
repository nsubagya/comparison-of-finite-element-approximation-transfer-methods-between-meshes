function u1=Lagrangeint(x,z,u,n)
for j=1:n
    for i=1:length(x)-1
        if z(j)>=x(i) && z(j)<x(i+1)
            I(j)=i;%element y1 in and Ix1=[x(I1) x(I1+1)], u1=[u(I1) u(I1+1)]
        end
    end
end
%boundary nodes
a=x(I);
b=x(I+1);

lambda=(z-a)./(b-a);
%nodal values of gauss points
u1=lambda.*u(I+1)+(1-lambda).*u(I);