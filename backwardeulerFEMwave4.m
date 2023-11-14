%solving du/dt+alpha*d u/dt=0 with Dirichlet boundry conditions using
%backward Euler and FE method with SUPG correction
%N-total time
%deltat-time step size
%x-FE space 1
%y-FE space 2
%swch- switch between L2projection (0) and lagrange interpolation (1) or
%Mass lumping (2)
%dtau- SUPG rate
%n-simpconser=0, simpconser+L2=1, simpconser+Masslump=2 (when swch=3)
%as an example use this "backwardeulerFEMwave4(0:0.1:10,0:0.1:10,0.1,1,0)"
%for non overlapping vertices use geny.m
function v2=backwardeulerFEMwave4(x,y,dt,tend,dtau,swch,n,gausspt)
%initial un in space x (t=0)
% u1=zeros(1,length(x));
% m=round(length(x)/2);
% u1(m-1)=1; u1(m)=1; u1(m+1)=1;
u1=uexact(x,0);
%save x to plot exact solution.
xinitial=x;
xexact=x(1):0.0001:x(length(x));
%variable for SUPG (dtau=0 means FE method without SUPG and dtau=1 means SUPG stabilization)
%dtau=0;
%piclate number
%picnum=8.84;
%boundary conditions
t=dt;
%ufirst=uexact(0,0);
%ufirst=0;
%ulast=0;
% hold on
% plot(x,u1,'k')


while t<=tend
    
    %initiate the mass matrix and stiffness matrix
    M=sparse(1:length(y),1:length(y),0,length(y),length(y));
    K=sparse(1:length(y),1:length(y),0,length(y),length(y));
    B=sparse(1:length(y),1:length(y),0,length(y),length(y));
    b=zeros(1,length(y));
    alphaq1=alpha(t);
    %finding mass and stiff matrices. this changes in every time step
    for i=1:length(y)-1
        H=y(i+1)-y(i);
        %alphaq1=(alpha(i)+alpha(i+1))/2;
        tau=dtau*H/alphaq1;
        %tau=(H/(2*alphaq1))*(1/tanh(picnum)-1/picnum)*dtau;
        M(i:i+1,i:i+1)=M(i:i+1,i:i+1)+(H/6)*[2 1;1 2];
        K(i:i+1,i:i+1)=K(i:i+1,i:i+1)+(1/2)*[0 -1;1 0];
        B(i:i+1,i:i+1)=B(i:i+1,i:i+1)+(1/H)*[1 -1;-1 1];
        %tauK=tau*K;
        tauB=tau*B;
    end
    %DONE with LHS
    A=M+dt*(alphaq1*K+(alphaq1^2)*tauB);
    LargeVal=1.0e8;
    %Dirichlet BC (we only need one boundry condition)
    A(1,1)=LargeVal;
    %A(length(A),length(A))=LargeVal;
    
    
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
    b=M*v1';%+tauK*v1';
    %we only need one boundry condition
    b(1)=LargeVal*uexact(0,t);
    %b(length(b))=LargeVal*ulast;
    %solve for new vn+1 in yn+1
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
%plot exact solution (%here t=tend when ploting)
%uexactplot=uexact(x,tend);
%uexactplot=uexact(xinitial,tend);
uinitialplot=uexact(xinitial,tend);
uexactplot=uexact(xexact,tend);
hold on
%plot(xinitial,uexactplot,'r',x,v2,'b')
plot(xexact,uexactplot,'r',xinitial,uinitialplot,'k',x,v2,'b')
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