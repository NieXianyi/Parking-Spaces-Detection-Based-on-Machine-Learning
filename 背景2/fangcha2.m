clear;
close all;
%% 用户定义
H = 2; %摄像机高度
L = 4; %摄像机深度
l = 2;  %车的宽度
h = 1.5; %车的高度


up=[195,435;
    338,435;
    481,435;
   618,435];
down=[1,542;
    278,540;
    548,538;
    805,536];



threshold = 20; %车辆的面积占车位面积的百分比的阈值
%threshold2 = 15; %车辆的面积占车位面积的百分比的阈值
%% 读取图形以及预处理
%A = (imread('空车位(扩展)(亮度调整).bmp'));
A = (imread('背景.jpg'));
A_2 = (imread('6.jpg'));
A = rgb2gray(A);
A_1_1 = im2bw(A,0.7);
A_2_2 = im2bw(rgb2gray(A_2),0.7);       %有车图的二值化图
%A_minus = A_2_2 - A_1_1;
A_minus = abs(rgb2gray(A_2) - A);            %灰度值差
A_minus2 = A_minus;
%thresh = graythresh(A_minus);
A_minus = im2bw(A_minus,0);

subplot(221);
imshow(A);
title('停车场');
subplot(222);
imshow(A_2);


title('已停放车辆的停车场');
subplot(2,2,4);
imshow(A_minus2);
title('差影运算');
%% 摄像机标定
theta1 = atan((H-h)/L);
demarcate_l = sin(theta1)*h/l; 
%得到视角在世界坐标系的比例

compensation_x = 0;         %一般不需要标定x方向
compensation_y = -round(demarcate_l*150);  %需要在Y方向补偿若干个像素
%% 计算出各个车位边线
for p = 1:1:size(down,1)
    for i = 1:1:down(p,2)-up(p,2)+1
        k = (up(p,1)-down(p,1))/(up(p,2)-down(p,2));
        B{p}(i,1) = round(up(p,1)+i*k);
        B{p}(i,2) = up(p,2)+i-1;
    end
end
    
%% 在原图标记出车位
 subplot(223);
 imshow(A);
 title('标记停车位');
 hold on
 for k = 1:length(B)
     %plot(B{k}(:,1), B{k}(:,2), 'G', 'LineWidth', 2)
 end
 for k = 1:length(B)-1
     %plot([B{k}(1,1),B{k+1}(1,1)],[B{k}(1,2),B{k+1}(1,2)], 'G', 'LineWidth', 2);
    % plot([B{k}(size(B{k},1),1),B{k+1}(size(B{k+1},1),1)],[B{k}(size(B{k},1),2),B{k+1}(size(B{k+1},1),2)], 'G', 'LineWidth', 2);
 end
%% 显示差影运算
figure;
subplot(221);
imshow(A_minus);
title('差影运算');
%% 闭运算
A_minus = bwmorph(A_minus,'close');% 闭运算
%% 显示闭运算
subplot(2,2,2);
imshow(A_minus);
title('闭运算');
% %% 在闭运算图标记停车位
subplot(2,2,3);
imshow(A_minus);
title('停车位标记');
hold on
for k = 1:length(B)
   % plot(B{k}(:,1), B{k}(:,2), 'G', 'LineWidth', 2)
end
for k = 1:length(B)-1
   % plot([B{k}(1,1),B{k+1}(1,1)],[B{k}(1,2),B{k+1}(1,2)], 'G', 'LineWidth', 2);
  %  plot([B{k}(size(B{k},1),1),B{k+1}(size(B{k+1},1),1)],[B{k}(size(B{k},1),2),B{k+1}(size(B{k+1},1),2)], 'G', 'LineWidth', 2);
end
%% 扫描每个像素检测是否有车
subplot(2,2,4);
imshow(A_2);
title('标记结果');%近似认为遮挡部分与未遮挡部分大小相同
hold on
SUM_BG = 0;%SUM_CAR = 0;%背景与车像素点
 car_judge = [];
for p = 1:1:length(B)-1          %看有没有车
    B1 = max(B{p}(:,2)) - min(B{p}(:,2));
    B2 = max(B{p+1}(:,2)) - min(B{p+1}(:,2));   %寻找最大最小值之差
    for i = 1:1:min(B1,B2)        %从一循环到最小的值
        for j = B{p}(i,1):1:B{p+1}(i,1)     %左边线遍历到右边
            SUM_BG = SUM_BG + 1;            %像素个数
            car_judge(SUM_BG) = A_minus2(B{p}(1,2)+i-1,j);
           % plot([j,j+1],[B{p}(1,2)+i-1,B{p}(1,2)+i-1], 'B', 'LineWidth', 2);
            if A_minus(B{p}(1,2)+i-1,j)== 1   %识别
%                 SUM_CAR = SUM_CAR + 1;
                %plot([j,j+1],[B{p}(1,2)+i-1,B{p}(1,2)+i-1], 'R', 'LineWidth', 2);                
            end
        end
    end
    temp = (std(car_judge(1:SUM_BG))).^2;
    if temp > threshold
        disp(['未标定时第',num2str(p),'个停车位：有车']);
    else
        disp(['未标定时第',num2str(p),'个停车位：无车']);
    end
    disp(['未标定时方差为',num2str(temp)]);
    SUM_BG = 0;
end

%% 分割线

%——————————————————————————————————————————%

%% 显示标定的停车位识别区域
figure
 subplot(221);
 imshow(A);
title('标定的停车位标记(检测区域)');
hold on
for k = 1:length(B)
   % plot(B{k}(:,1), B{k}(:,2), 'G', 'LineWidth', 2)
end
for k = 1:length(B)-1
   % plot([B{k}(1,1),B{k+1}(1,1)],[B{k}(1,2),B{k+1}(1,2)], 'G', 'LineWidth', 2);
   % plot([B{k}(size(B{k},1),1),B{k+1}(size(B{k+1},1),1)],[B{k}(size(B{k},1),2),B{k+1}(size(B{k+1},1),2)], 'G', 'LineWidth', 2);
end
for k = 1:length(B)
   % plot(B{k}(:,1)+compensation_x, B{k}(:,2)+compensation_y, 'Color',[1,0,1], 'LineWidth', 2)
end
for k = 1:length(B)-1
    %plot([B{k}(1,1)+compensation_x,B{k+1}(1,1)+compensation_x],[B{k}(1,2)+compensation_y,B{k+1}(1,2)+compensation_y], 'Color',[1,0,1], 'LineWidth', 2);
    %plot([B{k}(size(B{k},1)+compensation_x,1),B{k+1}(size(B{k+1},1)+compensation_x,1)],[B{k}(size(B{k},1)+compensation_y,2),B{k+1}(size(B{k+1},1)+compensation_y,2)], 'Color',[1,0,1], 'LineWidth', 2);
end
% 显示标定后的停车位区域
subplot(222);
imshow(A);
title('标定后的停车位标记(检测区域)');
hold on
for k = 1:length(B)
    %plot(B{k}(:,1)+compensation_x, B{k}(:,2)+compensation_y, 'G', 'LineWidth', 2)
end
for k = 1:length(B)-1
    %plot([B{k}(1,1)+compensation_x,B{k+1}(1,1)+compensation_x],[B{k}(1,2)+compensation_y,B{k+1}(1,2)+compensation_y], 'G',[1,0,1], 2);
   % plot([B{k}(size(B{k},1)+compensation_x,1),B{k+1}(size(B{k+1},1)+compensation_x,1)],[B{k}(size(B{k},1)+compensation_y,2),B{k+1}(size(B{k+1},1)+compensation_y,2)], 'G',[1,0,1], 2);
end
% 显示标定后的显示结果
subplot(2,2,3.5);
imshow(A_2);
title('标定后标记结果');%近似认为遮挡部分与未遮挡部分大小相同
hold on
SUM_BG = 0;SUM_CAR = 0;%背景与车像素点

for p = 1:1:length(B)-1          %看有没有车
    B1 = max(B{p}(:,2)) - min(B{p}(:,2));
    B2 = max(B{p+1}(:,2)) - min(B{p+1}(:,2));   %寻找最大最小值之差
    for i = 1:1:min(B1,B2)        %从一循环到最小的值
        for j = B{p}(i,1):1:B{p+1}(i,1)     %左边线遍历到右边
            SUM_BG = SUM_BG + 1;
            car_judge(SUM_BG) = A_minus2(B{p}(1,2)+i-1+compensation_y,j+compensation_x);
            %plot([j+compensation_x,j+1+compensation_x],[B{p}(1,2)+i-1+compensation_y,B{p}(1,2)+i-1+compensation_y], 'B', 'LineWidth', 2);
            if A_minus(B{p}(1,2)+i-1+compensation_y,j+compensation_x) == 1   %加入补偿
                SUM_CAR = SUM_CAR + 1;
              %  plot([j+compensation_x,j+1+compensation_x],[B{p}(1,2)+i-1+compensation_y,B{p}(1,2)+i-1+compensation_y], 'R', 'LineWidth', 2);                
            end
        end
    end
    temp = (std(car_judge(1:SUM_BG))).^2;
    if temp > threshold
        disp(['第',num2str(p),'个停车位：有车']);
    else
        disp(['第',num2str(p),'个停车位：无车']);
    end
    disp(['标定后方差为',num2str(temp)]);
    SUM_BG = 0;
    SUM_CAR = 0;
end