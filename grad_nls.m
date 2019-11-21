function [g,g1] = grad_nls(X,x,r,sigma2)
L = size(X,2);
t1 = 0; t2 = 0;
t11 = 0; t22 = 0;
ds = sum((x*ones(1,L)-X).^2,1);
ds = ds';
for i=1:L
    t1 = t1 + (1/sigma2(i))*((r(i)-ds(i)^(0.5))*(x(1)-X(1,i))/ds(i)^(0.5));
    t2 = t2 + (1/sigma2(i))*((r(i)-ds(i)^(0.5))*(x(2)-X(2,i))/ds(i)^(0.5));
    
    t11 = t11 + ((r(i)-ds(i)^(0.5))*(x(1)-X(1,i))/ds(i)^(0.5));
    t22 = t22 + ((r(i)-ds(i)^(0.5))*(x(2)-X(2,i))/ds(i)^(0.5));
end
g = -2*[t1; t2];
g1 = -2*[t11; t22];
end