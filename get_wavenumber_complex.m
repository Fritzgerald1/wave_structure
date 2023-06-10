function [RealK,ImagK,Omega,t] = get_wavenumber_complex(w_sca,Fun)
%%
% 输入量：
%	w_sca为无量纲数
%		w_sca = w*h/CT
% 输出量：
%	RealK也是无量纲数
%		RealK = k*h
%	t是模态数量
%	Omega是t个w_sca

%% 变量初始化
RealK = [];
ImagK = [];
Omega = [];
t=0;

%% 搜索波数解的范围
Ary = -1e-2;
Bry = 7;
dry = 1e-2;

Aiy = -1e-2;
Biy = 7;
diy = 1e-2;

%% 搜索波数解
parfor ii = 1:numel(w_sca) % 扫描频率的区间
    aw = w_sca(ii);
    for  ary = Ary:10*dry:Bry % 搜索波数解实部的范围
        for aiy = Aiy:10*diy:Biy
            s1=inf;
            sw1=0;
            sry1=0;
            siy1=0;
            bry = ary+11*dry; % 搜索范围的实部右边界
            biy = aiy+11*diy; % 搜索范围的虚部右边界
            %% 在小网格里，找到极小值解(syr1,siy1,sw1)
            for ry1 = ary:dry:bry
                for iy1 = aiy:diy:biy
                    y1 = ry1+1i*iy1;
                    h1 = Fun(y1,aw); % 频散方程行列式
                    hh1 = abs(h1);
                    if hh1 < s1
                        s1 = hh1;
                        sw1 = aw;
                        sry1 = ry1;
                        siy1 = iy1;
                    end
                end
            end

            if (abs(sry1-ary)>0.1*dry) && (abs(sry1-bry)>0.1*dry) && (abs(siy1-biy)>0.1*diy) && (abs(siy1-aiy)>0.1*diy) % 不在外边界上
                %% 当极小值不在边界时，判断该点是否是零点
                [zry,ziy,zero_flag]=Determine_zero_point_complex(sry1,siy1,sw1,dry,diy,Fun); % 判断零点，ssy1为每次迭代的解
                if zero_flag == 1 % 该点是零点
                    Omega=[Omega sw1];
                    RealK=[RealK zry];
                    ImagK=[ImagK ziy];
                    t=t+1;
                end
            end
        end
    end
end