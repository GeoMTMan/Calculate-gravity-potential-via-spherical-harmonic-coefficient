%% 本程序用于绘制重力势图
%  
%  李泯 编写于 2018/11/10
% 
%% 程序初始化
clear; 
close all;
clc;

%% 导入数据
load CSnm; % 重力势系数
load coast; % 海岸线数据

%% 参数定义与初始化
nmax = 120; % 纬度方向网格数，经度方向网格数为 2*nmax
order_use = 12; % 要计算的阶数
GM = 3.986004415*10^14; % m^2/sec^3
R = 6378136.3; % m
% 参数初始化
theta = linspace(0, pi, nmax+1); % 纬度方向，球坐标系坐标
phi = linspace(0, 2*pi, 2*nmax+1); % 经度方向，球坐标系坐标
[long_v,lat_v]=meshgrid(phi*180/pi,90-theta*180/pi); % 生成经纬度网格，平面坐标
V = zeros(length(theta), length(phi)); % 重力势

%% 计算重力势
V = V + Cnm(1);
for n = 1:order_use
    Pnm = legendre(n, cos(theta),'norm');
    for m = 0:n
        V = V + Pnm(m+1,:)'*(Cnm((1+n)*n/2+m+1)*cos(m*phi) + ...
            Snm((1+n)*n/2+m+1)*sin(m*phi));
    end
end
V = -(GM/R)*V;
%% 绘图
% 绘制平面图
figure; worldmap world;
pcolorm(lat_v,long_v,V);
plotm(lat,long,'w-');
colorbar('location','EastOutside');
% 绘制三维图
figure; sphere;
h = findobj(gcf, 'Type', 'surface');
set(h, 'CData', V, 'FaceColor', 'texturemap');
axis equal; axis off;
colorbar('location','EastOutside');
% 旋转三维视图
%for i=0:360
%    view(i,30);
%    drawnow;
%end
