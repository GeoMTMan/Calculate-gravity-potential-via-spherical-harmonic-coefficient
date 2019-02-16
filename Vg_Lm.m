%% ���������ڻ���������ͼ
%  
%  ���� ��д�� 2018/11/10
% 
%% �����ʼ��
clear; 
close all;
clc;

%% ��������
load CSnm; % ������ϵ��
load coast; % ����������

%% �����������ʼ��
nmax = 120; % γ�ȷ��������������ȷ���������Ϊ 2*nmax
order_use = 12; % Ҫ����Ľ���
GM = 3.986004415*10^14; % m^2/sec^3
R = 6378136.3; % m
% ������ʼ��
theta = linspace(0, pi, nmax+1); % γ�ȷ���������ϵ����
phi = linspace(0, 2*pi, 2*nmax+1); % ���ȷ���������ϵ����
[long_v,lat_v]=meshgrid(phi*180/pi,90-theta*180/pi); % ���ɾ�γ������ƽ������
V = zeros(length(theta), length(phi)); % ������

%% ����������
V = V + Cnm(1);
for n = 1:order_use
    Pnm = legendre(n, cos(theta),'norm');
    for m = 0:n
        V = V + Pnm(m+1,:)'*(Cnm((1+n)*n/2+m+1)*cos(m*phi) + ...
            Snm((1+n)*n/2+m+1)*sin(m*phi));
    end
end
V = -(GM/R)*V;
%% ��ͼ
% ����ƽ��ͼ
figure; worldmap world;
pcolorm(lat_v,long_v,V);
plotm(lat,long,'w-');
colorbar('location','EastOutside');
% ������άͼ
figure; sphere;
h = findobj(gcf, 'Type', 'surface');
set(h, 'CData', V, 'FaceColor', 'texturemap');
axis equal; axis off;
colorbar('location','EastOutside');
% ��ת��ά��ͼ
%for i=0:360
%    view(i,30);
%    drawnow;
%end
