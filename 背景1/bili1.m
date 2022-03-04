%% 第一组

clc;
close all;
%% 用户定义
H = 15; %摄像机高度
L = 30; %摄像机深度
l = 2;  %车的宽度
h = 1.5; %车的高度

up = [210,470;          %上面若干点
    352,469;
    502,468;
    652,466];
down = [12,588;        %下面若干点
    280,589;
    560,589;
    830,590];

%threshold = 20; %车辆的面积占车位面积的百分比的阈值
%threshold2 = 15; %车辆的面积占车位面积的百分比的阈值

%% 摄像机标定
theta1 = atan((H-h)/L);
demarcate_l = sin(theta1)*h/l; 
%得到视角在世界坐标系的比例

compensation_x = 0;         %一般不需要标定x方向
compensation_y = round(demarcate_l*150);  %需要在Y方向补偿若干个像素

%% 读取图形以及预处理
%A = (imread('空车位(扩展)(亮度调整).bmp'));
A = (imread('没有车.jpg'));
A = rgb2gray(A);
A_1 = (imread('12.jpg'));
A_1 = rgb2gray(A_1);
A_2=A_1-A;
% A = rgb2gray(A);
% A_1_1 = im2bw(A,0.7);
% A_2_2 = im2bw(rgb2gray(A_2),0.7);       %有车图的二值化图
% %A_minus = A_2_2 - A_1_1;
% A_minus = abs(rgb2gray(A_2) - A);            %灰度值差
% A_minus2 = A_minus;
% %thresh = graythresh(A_minus);
% A_minus = im2bw(A_minus,0);


%% 显示差影运算
figure;

imshow(A_2);
title('差影运算');
[m,n]=size(A_2);

w=fspecial('gaussian',[5 5],1.6);
A_2=imfilter(A_2,w,'replicate');

z=0;
for i = (470+compensation_x):(590+compensation_y)
    for j = (12+compensation_x):(280)
        if A_2(i,j)>0
            z=z+1;
        end
    end
end
q=double(z/(m*n));

 disp(['第一个车位的比例为',num2str(q)]);
 
 z=0;
for i = (469+compensation_x):(589+compensation_y)
    for j = (352+compensation_x):(502)
        if A_2(i,j)>0
            z=z+1;
        end
    end
end
q=double(z/(m*n));
 disp(['第二个车位的比例为',num2str(q)]);
 
  z=0;
for i = (468+compensation_x):(589+compensation_y)
    for j = (560+compensation_x):(830)
        if A_2(i,j)>0
            z=z+1;
        end
    end
end
q=double(z/(m*n));
 disp(['第三个车位的比例为',num2str(q)]);