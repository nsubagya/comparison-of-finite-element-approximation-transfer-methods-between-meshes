%x-FE space 1
%y-FE space 2
%u-sulution set in x
%graphflag- if graphflag >0 graph is on O/W off
function u2=LagInt7(x,y,u)
%Example command to type in command window
%FEMparrun3([-0.5:0.03:0.5],[-0.5:0.09:-0.2,-0.2:0.05:0.2,0.2:0.09:0.5],5*rand(length([-0.5:0.03:0.5]),1),1)
z=union(x,y);
%initiating the output vector
u2=zeros(1,length(y));
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
    %interpolate only if z(k) point is in y
    n=ismember(z(k),y);
    switch n
        case 1
            %I=[i j];(keeping this to troubleshooting if needed)
            alphax=(y(j)-x(i))/(x(i+1)-x(i));
            u2(j)=alphax*u(i+1)+(1-alphax)*u(i);
    end
end
u2(length(y))=u(length(x));

end
