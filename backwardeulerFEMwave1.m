%solving du/dt+alpha*d u/dt=0 with Dirichlet boundry conditions using
%backward Euler and FE method
%N-total time
%deltat-time step size
%x-FE space 1
%y-FE space 2
%swch- switch between L2projection (0) and lagrange interpolation (1) or
%Mass lumping (2)
%as an example use this "backwardeulerFEMwave1(0:0.5:10,0:0.25:10,5,1)"
%for non overlapping vertices use geny.m
function v2=backwardeulerFEMwave1(x,y,dt,tend,swch)
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
%plot initial condition
% hold on
% plot(x,u1,'k')


while t<=tend
    
    %initiate the mass matrix and stiffness matrix
    M=sparse(1:length(y),1:length(y),0,length(y),length(y));
    K=sparse(1:length(y),1:length(y),0,length(y),length(y));
    b=zeros(1,length(y));
    alphaq1=alpha(t);
    %finding mass and stiff matrices. this changes in every time step
    for i=1:length(y)-1
        H=y(i+1)-y(i);
        %alphaq1=(alpha(i)+alpha(i+1))/2;
        M(i:i+1,i:i+1)=M(i:i+1,i:i+1)+(H/6)*[2 1;1 2];
        K(i:i+1,i:i+1)=K(i:i+1,i:i+1)+(1/2)*[0 -1;1 0];
        
    end
    %DONE with LHS
    A=M+dt*alphaq1*K;
    LargeVal=1.0e8;
    %Dirichlet BC (we only need one boundry condition)
    A(1,1)=LargeVal;
    %A(length(A),length(A))=LargeVal;
    
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
    b=M*v1';%+deltat*F';
    %we only need one boundry condition
    b(1)=LargeVal*uexact(0,t);
    % b(length(b))=LargeVal*ulast;
    %solve for new vn+1 in yn+1 (solve xA=b)
    v2=b'/A;
    
    
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
%plot exact solution (for final t=tend)
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