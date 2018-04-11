function [y] = mycummin(x)
minSoFar = x(1);
y=zeros(size(x));
y(1)=minSoFar;
for i=2:length(x)
    if(minSoFar>x(i))
        minSoFar = x(i);  
    end
    y(i)=minSoFar;    
end