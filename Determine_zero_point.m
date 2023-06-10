function  [true_sy,zero_flag]=Determine_zero_point(y,w,dy,Fun)
%% Determine_zero_point
% 判断(y,w)是否为零点
% 输出：
%	ssy1:末项为零点y的精确解；
%	DetF:每次迭代的行列式，当DetF(1)/DetF(end) > 1e7时
%	zero_flag:为1时，该点为零点；为0时，该点不是零点。

true_sy = []; % y的精确解
zero_flag = 0;
Detmin=inf;
sy = 0;
b = 20;% 迭代次数
step = 2;% 将子区间[yy-ty,yy+ty]划分为2*{step}份
DetF=zeros(1,b);% 存储每个迭代的行列式
ssy1=zeros(1,b);% 存储每个迭代的解

%% 搜索区间为[y-ty,y+ty]，步长(ty/step)
ay = y-dy;
by = y+dy;
dy = dy/step;

for k = 1:b
	for y1=ay:dy:by
		[h1,flag] = Fun(y1,w);
		hh1=abs(h1);
		if hh1<Detmin
			Detmin=hh1;
			sy=y1;
		end
	end
	DetF(k)=Detmin; % 存储每个迭代的行列式
	ssy1(k) = sy;% 存储每个迭代的解

	%% 缩小搜索区间
	ay = sy-dy;
	by = sy+dy;
	dy = dy/step;
end

if DetF(1)/DetF(end) > 1e2 % 为零点
	if flag
	zero_flag = 1;
	true_sy = ssy1(end);
end
end