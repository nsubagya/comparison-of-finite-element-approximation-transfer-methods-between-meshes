function u3=gauss(x,y,j,u,gausspt)
%initiate u2 (left and right mass of an element j)
u2=zeros(2,1);
%boundary of j element
a=y(j);
b=y(j+1);
%for changing to a,b from -1,1
dy=(b-a)/2;
%-1,1 nodes and waights
nfull=[0 0 0 0 0;-1/sqrt(3) 1/sqrt(3) 0 0 0;-sqrt(3/5) 0 sqrt(3/5) 0 0;-sqrt(3/7+2/7*sqrt(6/5)) -sqrt(3/7-2/7*sqrt(6/5)) sqrt(3/7-2/7*sqrt(6/5)) sqrt(3/7+2/7*sqrt(6/5)) 0;-1/3*sqrt(5+2*sqrt(10/7)) -1/3*sqrt(5-2*sqrt(10/7)) 0 1/3*sqrt(5-2*sqrt(10/7)) 1/3*sqrt(5+2*sqrt(10/7))];
wfull=[2 0 0 0 0;1 1 0 0 0;5/9 8/9 5/9 0 0;(18-sqrt(30))/36 (18+sqrt(30))/36 (18+sqrt(30))/36 (18-sqrt(30))/36 0;(322-13*sqrt(70))/900 (322+13*sqrt(70))/900 128/225 (322+13*sqrt(70))/900 (322-13*sqrt(70))/900];
% n=[-sqrt(3/5) 0 sqrt(3/5)];
% w=[5/9 8/9 5/9];
n=nfull(gausspt,1:gausspt);
w=wfull(gausspt,1:gausspt);
m=(b-a)/2*n+(a+b)/2;
%u,phi values at points n
ug=Lagrangeint(x,m,u,length(n));
phigright=Lagrangeint([y(j) y(j+1)],m,[0 1],length(n));
phigleft=Lagrangeint([y(j) y(j+1)],m,[1 0],length(n));
for i=1:length(n)
   u2(1)=u2(1)+w(i)*ug(i)*phigleft(i);
   u2(2)=u2(2)+w(i)*ug(i)*phigright(i);
end
%mass of nodes
u3=dy*u2;