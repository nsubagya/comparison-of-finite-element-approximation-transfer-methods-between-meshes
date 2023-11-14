%x-FE space 1
%y-FE space 2
%u-sulution set in x (pld mesh)
%U-continues nodal values of the new mesh
%n-simpconser=0, simpconser+L2=1, simpconser+Masslump=2
%swch-Alauzet=0, L2 proj=1, Lagint=2, masslump=3
%mend-number of times to run the script
function U=runningscript(x,y,u,n,mend,swch,gausspt)
%Example command to type in command window
%U=runningscript(0:0.5:10,0:0.25:10,exp(-(x-5).^2),0,1000,0)
dm=1;
m=dm;
tic
while m<=mend
    switch swch
        case 0
            U=Alauzet(x,y,u,n);
        case 1
            U=L2proj7(x,y,u);
        case 2
            U=LagInt7(x,y,u);
        case 3
            U=MassLump7(x,y,u);
        case 4
            U=L2projwosupmesh8(x,y,u,gausspt);
    end
    xold=x;
    x=y;
    y=xold;
    u=U;
    m=m+dm;
end
toc
end