function [error,varargout] = lamb_sym_real2(kd,wd)
%% 材料属性
CL = 6.35;
CT = 3.13;
h=0.5;

%% 量纲还原
k = kd/h;
w = wd*CT/h;
p = sqrt(w.^2/CL^2 - k.^2);
q = sqrt(w.^2/CT^2 - k.^2);

%% 频散方程
M(1,1) = 2i*k*p*sin(p*h);
M(1,2) = (k^2-q^2)*sin(q*h);
M(2,1) = (q^2-k^2)*cos(p*h);
M(2,2) = -2i*k*q*cos(q*h);
error = det(M);

if nargout == 2
	flag = 1;
	if abs(p)<1e-1
		flag = 0;
	end
	if abs(q)<1e-1
		flag =0;
	end
	varargout{1} = flag;
end
end
