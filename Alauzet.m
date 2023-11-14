%x-FE space 1
%y-FE space 2
%u-sulution set in x (pld mesh)
%p1,p2-element wise approxiation with first and last nodal value of the
%element
%U-continues nodal values of the new mesh
%n-simpconser=0, simpconser+L2=1, simpconser+Masslump=2
function U=Alauzet(x,y,u,n)
%Example command to type in command window
%U=Alauzet(0:0.5:10,0:0.25:10,exp(-(x-5).^2))
%construct the super mesh
z=union(x,y);
%initiating the output vector gradiant and mass
u3=zeros(1,length(y)-1);
gradunew=zeros(1,length(y)-1);
b=zeros(1,length(y));
M=sparse(1:length(y),1:length(y),0,length(y),length(y));
Mlump=sparse(1:length(y),1:length(y),0,length(y),length(y));%maybe Mlump=M but not sure
umax=ones(length(y)-1,1)*-1.0e8;
umin=ones(length(y)-1,1)*1.0e8;
%mass matrix and lumped mass matrix
for i=1:length(y)-1
    H=y(i+1)-y(i);
    M(i:i+1,i:i+1)=M(i:i+1,i:i+1)+(H/6)*[2 1;1 2];
    Mlump(i:i+1,i:i+1)=Mlump(i:i+1,i:i+1)+(H/2)*[1 0;0 1];
end

%find interval values (starting value) of the y interval
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
    %finding maximum and minimum value of mesh i wrt j element
    umax(j)=max([u(i),u(i+1),umax(j)]);
    umin(j)=min([u(i),u(i+1),umin(j)]);
    %finding mesh element size
    hx=x(i+1)-x(i);
    hy(j)=y(j+1)-y(j); %in some cases this counts twice
    hz=z(k+1)-z(k);
    %area divition factor %area of current part of the y element %devide to find hight of the element A=h*l
    u2=hz*(u(i)+u(i+1))/(2*hy(j));
    %A1/hy + A2/hy =(A1+A2)/hy hight of the total y element (take this as the mass at the bericenter)
    u3(j)=u3(j)+u2;
    %finding the gradiant of the element in x mesh (since linear it will same in z and y values in x)
    gradu=(u(i+1)-u(i))*hz/(hy(j)*hx);
    %A1/hy + A2/hy =(A1+A2)/hy gradiant change of the y element
    gradunew(j)=gradunew(j)+gradu;
    
    
    %end
end
%find first and last value of the element
%to use Alauzet's formula find first element and the node
for j=1:length(y)-1
    p1(j)=u3(j)+gradunew(j)*-hy(j)/2;
    p2(j)=u3(j)+gradunew(j)*hy(j)/2;
    
    up1=min(p1(j),p2(j));
    up2=max(p1(j),p2(j)); %maybe we dont need to find both
    
    %for p1 is the minimum
    if up1==p1(j)
        pM2=min(up2,umax(j));
        pM1=up1+(up2-pM2);
        
        p1(j)=max(pM1,umin(j));
        p2(j)=pM2-(p1(j)-pM1);
    %for p2 is the minimum
    else
        pM1=min(up2,umax(j));
        pM2=up1+(up2-pM1);
        
        p2(j)=max(pM2,umin(j));
        p1(j)=pM1-(p2(j)-pM2);
    end
    switch n
        case 0
            switch j
                case 1
                    U(j)=(hy(j)*p1(j))/(hy(j));
                    %find the other nodes (vertices) except last one
                otherwise
                    U(j)=(hy(j-1)*p2(j-1)+hy(j)*p1(j))/(hy(j)+hy(j-1));
            end
        otherwise
            b(j)=b(j)+(p1(j)+4*(p1(j)+p2(j))/4)*(y(j+1)-y(j))/6;
            b(j+1)=b(j+1)+(p2(j)+4*(p1(j)+p2(j))/4)*(y(j+1)-y(j))/6;
    end
end
%last nodel value
switch n
    case 0
        U(length(y))=(hy(length(y)-1)*p2(length(y)-1))/(hy(length(y)-1));
    case 1
        U=b/M;
    case 2
        U=b/Mlump;
end
end
