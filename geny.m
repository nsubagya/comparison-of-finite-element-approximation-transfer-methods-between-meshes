%to get [0 0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10] use w= 0:0.5:10 and
%x=0:10.
function y=geny(w,x)
if x(1)~=w(1)
    fprintf('x, w end pont values must be same (check x(1)=w(1)?)')
elseif x(length(x))~=w(length(w))
    fprintf('x, w end pont values must be same (check x(last)=w(last)?)')
else
B=intersect(w,x);
y=unique([x(1) setdiff(w,B) x(length(x))]);
end