%solving du/dt=d^2 u/dt^2 +f with Dirichlet boundry conditions (Heat equation)
%N-total time
%deltat-time step size
%x-FE space 1
%y-FE space 2
%swch- switch between L2projection (0) and lagrange interpolation (1) or
%Mass lumping (2) and Alauzeth (3) L2 projection without supermesh (4)
%n-simpconser=0, simpconser+L2=1, simpconser+Masslump=2 (when swch=3)
%as an example use this "backwardeulerFEM1(0:0.5:10,0:0.25:10,5,1)"
%for non overlapping vertices use geny.m
function v2=backwardeulerFEM1(x,y,dt,tend,swch,n,gausspt)
%initial un in space x (t=0)
%u1=5*(-((x/5)-1).^2 +1); %un in xn
u1=uexact(x,0);
xinitial=x;
xexact=x(1):0.0001:x(length(x));
% u1=zeros(1,length(x));
% m=round(length(x)/2);
% u1(m-1)=1; u1(m)=1; u1(m+1)=1;
t=dt;
%boundary conditions
% ufirst=0;
% ulast=0;
% hold on
% plot(xinitial,u1,'m')


while t<=tend
    
    %initiate the mass matrix and stiffness matrix
    M=sparse(1:length(y),1:length(y),0,length(y),length(y));
    K=sparse(1:length(y),1:length(y),0,length(y),length(y));
    b=zeros(1,length(y));
    F=zeros(1,length(y));
    %finding mass and stiff matrices. this changes in every time step
    for i=1:length(y)-1
        H=y(i+1)-y(i);
        alphaq1=(alpha(i)+alpha(i+1))/2;
        M(i:i+1,i:i+1)=M(i:i+1,i:i+1)+(H/6)*[2 1;1 2];
        K(i:i+1,i:i+1)=K(i:i+1,i:i+1)+(alphaq1/H)*[1 -1;-1 1];
        F(i:i+1)=F(i:i+1)+(H/2)*[f(i,t) f(i+1,t)];
        
    end
    %DONE with LHS
    A=M+dt*K;
    LargeVal=1.0e8;
    %Dirichlet BC
    A(1,1)=LargeVal;
    A(length(A),length(A))=LargeVal;
    
    %convert un in xn to vn in xn+1 (or say yn)
    %u1=proj(x,y,u1);
    switch swch
        case 0 %L2 projection
            v1=L2proj7(x,y,u1);
        case 1 %Lagrangian interpolation
            v1=LagInt7(x,y,u1);
        case 2 %Mass Lumping
            v1=MassLump7(x,y,u1);
        case 3 %Alauzet
            v1=Alauzet(x,y,u1,n);
        case 4
            v1=L2projwosupmesh8(x,y,u1,gausspt);
    end
    
    %find RHS
    b=M*v1'+dt*F';
    b(1)=LargeVal*uexact(0,t);
    b(length(b))=LargeVal*uexact(xinitial(length(xinitial)),t);
    %solve for new vn+1 in yn+1
    v2=b'/A;
    
    
    %plot the results(in every time step)
    %figure
    %     hold on
    %     plot(y,v2)
    %1 ST TIME STEP IS DONE!!!!!
    
    %setting values for next time step
    u1=v2;    %this will act as vn in yn
    %save xn old and yn old
    xold=x;
    
    %change the FE spaces for next time step
    x=y;
    y=xold;
    
    %     uexact=5*sin(pi*x/1).*exp(-n*(pi/1).^2);
    %     hold on
    %     plot(x,uexact,'k')
    t=t+dt;
end
%plot exact solution
uinitialplot=uexact(xinitial,tend);
uexactplot=uexact(xexact,tend);
hold on
%plot(xinitial,uexactplot,'r',x,v2,'b')
plot(xexact,uexactplot,'r',xinitial,uinitialplot,'k',x,v2,'b')
end
function uexact=uexact(x,t)
%uexact=exp(-((x-2)-alpha*t).^2);
uexact=5*sin(pi*x/10).*exp(-t*(pi/10).^2);
%uexact=((x-alpha*t)>4).*((x-alpha*t)<7);
end
function alpha=alpha(x)
alpha=1;
end
function f=f(x,t)
f=0;
end