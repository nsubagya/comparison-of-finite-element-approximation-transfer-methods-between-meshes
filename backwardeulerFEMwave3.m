%solving du/dt+alpha*d u/dt=0 with Dirichlet boundry conditions and using
%backward Euler and backward method in space
%N-total time
%deltat-time step size
%x-FE space 1
%y-FE space 2
%swch- switch between L2projection (0) and lagrange interpolation (1) or
%Mass lumping (2)
%as an example use this "backwardeulerFEMwave3(0:0.5:10,0:0.25:10,5,1)"
%for non overlapping vertices use geny.m
function v2=backwardeulerFEMwave3(x,y,dt,tend,swch)
%initial un in space x (t=0)
% u1=zeros(1,length(x));
% m=round(length(x)/2);
% u1(m-1)=1; u1(m)=1; u1(m+1)=1;
u1=uexact(x,0);
%save x to plot exact solution.
xinitial=x;
%time step starting value
t=dt;
%boundary conditions
ufirst=uexact(0,0);
%ufirst=0;
%ulast=0;
% hold on
% plot(x,u1,'r')


while t<=tend
    
    %initiate the mass matrix and stiffness matrix
    M=sparse(1:length(y),1:length(y),0,length(y),length(y));
    b=zeros(1,length(y));
    %finding mass and stiff matrices. this changes in every time step
    for i=2:length(y)
        H=y(i)-y(i-1);
        %alphaq1=(alpha(i)+alpha(i+1))/2;
        lambda=dt/H;
        M(i,i-1:i)=M(i,i-1:i)+[-lambda*alpha,1+lambda*alpha];
        
    end
    
    
    LargeVal=1.0e8;
    %Dirichlet BC (we only need one boundry condition)
    M(1,1)=LargeVal;
    %     M(length(M),length(M))=LargeVal;
    %M=M(2:length(M)-1,2:length(M)-1);
    %DONE with LHS
    
    
    %convert un in xn to vn in xn+1 (or say yn)
    %u1=proj(x,y,u1);
    switch swch
        case 0 %L2 projection
            v1=L2proj7(x,y,u1);
        case 1 %Lagrangian interpolation
            v1=LagInt7(x,y,u1);
        case 2 %Mass Lumping
            v1=MassLump7(x,y,u1);
    end
    
    %find RHS
    b=v1;%+deltat*F';
    %we only need one boundry condition
    b(1)=LargeVal*uexact(0,t);
    %b(length(b))=LargeVal*ulast;
    %b=b(2:length(b)-1);
    %solve for new vn+1 in yn+1
    v2=M\b';
    
    
    %plot the results(in every time step)
    %     figure
    %         hold on
    %         plot(y,v2)
    %1 ST TIME STEP IS DONE!!!!!
    
    %setting values for next time step
    u1=v2;    %this will act as vn in yn
    %save xn old and yn old
    xold=x;
    
    %change the FE spaces for next time step
    x=y;
    y=xold;
    %exact solution for every time step
    %     uexact=5*sin(pi*x/1).*exp(-n*(pi/1).^2);
    %     hold on
    %     plot(x,uexact,'k')
    t=t+dt;
end
%plot exact solution
uexactplot=uexact(xinitial,tend);
hold on
plot(xinitial,uexactplot,'r')

hold on
plot(x,v2)
end
function uexact=uexact(x,t)
uexact=exp(-((x-2)-alpha*t).^2);
%uexact=((x-alpha*t)>0.4).*((x-alpha*t)<1.6);
end
function alpha=alpha(t)
alpha=1;
end
function f=f(x,t)
f=0;
end