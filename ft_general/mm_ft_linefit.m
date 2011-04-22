function [m,b] = linefit(xdata, ydata, num, colstr)
% [m,b] = linefit(xdata, ydata, num, colstr)
%

degree = 1;
MB = polyfit(xdata,ydata,degree);
m = MB(1);
b = MB(2);
%fprintf('m%s = %.6f, b%s = %.6f\n',int2str(num),MB(1),int2str(num),MB(2))
X = [min(xdata) max(xdata)];
Ms = ones(1,2)*MB(1);
Bs = ones(1,2)*MB(2);
Y = Ms.*X + Bs;
plot(X,Y,colstr,'LineWidth',2);

end

