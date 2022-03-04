function varargout = jiemian(varargin)
% JIEMIAN MATLAB code for jiemian.fig
%      JIEMIAN, by itself, creates a new JIEMIAN or raises the existing
%      singleton*.
%
%      H = JIEMIAN returns the handle to a new JIEMIAN or the handle to
%      the existing singleton*.
%
%      JIEMIAN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in JIEMIAN.M with the given input arguments.
%
%      JIEMIAN('Property','Value',...) creates a new JIEMIAN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before jiemian_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to jiemian_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help jiemian

% Last Modified by GUIDE v2.5 18-Dec-2018 21:07:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @jiemian_OpeningFcn, ...
                   'gui_OutputFcn',  @jiemian_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before jiemian is made visible.
function jiemian_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to jiemian (see VARARGIN)

% Choose default command line output for jiemian
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes jiemian wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = jiemian_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in open.
function open_Callback(hObject, eventdata, handles)
% hObject    handle to open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
image = imread('2.jpg');
imshow(image);

% --- Executes on button press in shibie.
function shibie_Callback(hObject, eventdata, handles)
% hObject    handle to shibie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% 用户定义
H = 2; %摄像机高度
L = 4; %摄像机深度
l = 2;  %车的宽度
h = 1.5; %车的高度

up = [205,467;          %上面若干点
    350,470;
    505,470;
    655,470];
down = [10,590;        %下面若干点
    280,590;
    556,590;
    850,590];

threshold = 20; %车辆的面积占车位面积的百分比的阈值
%threshold2 = 15; %车辆的面积占车位面积的百分比的阈值
%% 读取图形以及预处理
%A = (imread('空车位(扩展)(亮度调整).bmp'));
A = (imread('没有车.jpg'));
A_2 = (imread('2.jpg'));
A = rgb2gray(A);
A_1_1 = im2bw(A,0.7);
A_2_2 = im2bw(rgb2gray(A_2),0.7);       %有车图的二值化图
%A_minus = A_2_2 - A_1_1;
A_minus = abs(rgb2gray(A_2) - A);            %灰度值差
A_minus2 = A_minus;
%thresh = graythresh(A_minus);
A_minus = im2bw(A_minus,0);

% subplot(221);
% imshow(A);
% title('停车场');
% subplot(222);
% imshow(A_2);


% title('已停放车辆的停车场');
% subplot(2,2,4);
% imshow(A_minus2);
% title('差影运算');
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
%  subplot(223);
%  imshow(A);
%  title('标记停车位');
 hold on
 for k = 1:length(B)
    plot(B{k}(:,1), B{k}(:,2), 'G', 'LineWidth', 2)
 end
 for k = 1:length(B)-1
     %plot([B{k}(1,1),B{k+1}(1,1)],[B{k}(1,2),B{k+1}(1,2)], 'G', 'LineWidth', 2);
     %plot([B{k}(size(B{k},1),1),B{k+1}(size(B{k+1},1),1)],[B{k}(size(B{k},1),2),B{k+1}(size(B{k+1},1),2)], 'G', 'LineWidth', 2);
 end
%% 显示差影运算
% figure;
% subplot(221);
% imshow(A_minus);
% title('差影运算');
%% 闭运算
A_minus = bwmorph(A_minus,'close');% 闭运算
%% 显示闭运算
 %subplot(2,2,2);
 %imshow(A_minus);
 %title('闭运算');
 %% 在闭运算图标记停车位
 %subplot(2,2,3);
 %imshow(A_minus);
 %title('停车位标记');
hold on
for k = 1:length(B)
   % plot(B{k}(:,1), B{k}(:,2), 'G', 'LineWidth', 2)
end
for k = 1:length(B)-1
    %plot([B{k}(1,1),B{k+1}(1,1)],[B{k}(1,2),B{k+1}(1,2)], 'G', 'LineWidth', 2);
    %plot([B{k}(size(B{k},1),1),B{k+1}(size(B{k+1},1),1)],[B{k}(size(B{k},1),2),B{k+1}(size(B{k+1},1),2)], 'G', 'LineWidth', 2);
end
%% 扫描每个像素检测是否有车
% subplot(2,2,4);
% imshow(A_2);
% title('标记结果');%近似认为遮挡部分与未遮挡部分大小相同
hold on
SUM_BG = 0;SUM_CAR = 0;%背景与车像素点
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
                 SUM_CAR = SUM_CAR + 1;
               % plot([j,j+1],[B{p}(1,2)+i-1,B{p}(1,2)+i-1], 'R', 'LineWidth', 2);                
            end
        end
    end
    temp = (std(car_judge(1:SUM_BG))).^2;
    SUM_BG = 0;
end

%% 显示标定的停车位识别区域
% figure
%  subplot(221);
%  imshow(A);
% title('标定的停车位标记(检测区域)');
hold on
for k = 1:length(B)
    %plot(B{k}(:,1), B{k}(:,2), 'G', 'LineWidth', 2)
end
for k = 1:length(B)-1
    plot([B{k}(1,1),B{k+1}(1,1)],[B{k}(1,2),B{k+1}(1,2)], 'G', 'LineWidth', 2);
    plot([B{k}(size(B{k},1),1),B{k+1}(size(B{k+1},1),1)],[B{k}(size(B{k},1),2),B{k+1}(size(B{k+1},1),2)], 'G', 'LineWidth', 2);
end
for k = 1:length(B)
    %plot(B{k}(:,1)+compensation_x, B{k}(:,2)+compensation_y, 'Color',[1,0,1], 'LineWidth', 2)
end
for k = 1:length(B)-1
    %plot([B{k}(1,1)+compensation_x,B{k+1}(1,1)+compensation_x],[B{k}(1,2)+compensation_y,B{k+1}(1,2)+compensation_y], 'Color',[1,0,1], 'LineWidth', 2);
    %plot([B{k}(size(B{k},1)+compensation_x,1),B{k+1}(size(B{k+1},1)+compensation_x,1)],[B{k}(size(B{k},1)+compensation_y,2),B{k+1}(size(B{k+1},1)+compensation_y,2)], 'Color',[1,0,1], 'LineWidth', 2);
end
% 显示标定后的停车位区域
% subplot(222);
% imshow(A);
% title('标定后的停车位标记(检测区域)');
hold on
for k = 1:length(B)
    plot(B{k}(:,1)+compensation_x, B{k}(:,2)+compensation_y, 'G', 'LineWidth', 2)
end
for k = 1:length(B)-1
    plot([B{k}(1,1)+compensation_x,B{k+1}(1,1)+compensation_x],[B{k}(1,2)+compensation_y,B{k+1}(1,2)+compensation_y], 'G',[1,0,1], 2);
    plot([B{k}(size(B{k},1)+compensation_x,1),B{k+1}(size(B{k+1},1)+compensation_x,1)],[B{k}(size(B{k},1)+compensation_y,2),B{k+1}(size(B{k+1},1)+compensation_y,2)], 'G',[1,0,1], 2);
end
% 显示标定后的显示结果
% subplot(2,2,3.5);
% imshow(A_2);
% title('标定后标记结果');%近似认为遮挡部分与未遮挡部分大小相同
hold on

SUM_BG = 0;SUM_CAR = 0;%背景与车像素点
%看有没有车

Z = [0,0,0];
for p = 1:1:length(B)-1          
    B1 = max(B{p}(:,2)) - min(B{p}(:,2));
    B2 = max(B{p+1}(:,2)) - min(B{p+1}(:,2));   %寻找最大最小值之差
    for i = 1:1:min(B1,B2)        %从一循环到最小的值
        for j = B{p}(i,1):1:B{p+1}(i,1)     %左边线遍历到右边
            SUM_BG = SUM_BG + 1;
            car_judge(SUM_BG) = A_minus2(B{p}(1,2)+i-1+compensation_y,j+compensation_x);
          % plot([j+compensation_x,j+1+compensation_x],[B{p}(1,2)+i-1+compensation_y,B{p}(1,2)+i-1+compensation_y], 'B', 'LineWidth', 2);
           
          if A_minus(B{p}(1,2)+i-1+compensation_y,j+compensation_x) == 1   %加入补偿
                SUM_CAR = SUM_CAR + 1;
             %plot([j+compensation_x,j+1+compensation_x],[B{p}(1,2)+i-1+compensation_y,B{p}(1,2)+i-1+compensation_y], 'R', 'LineWidth', 2);                
          end
            
        end
    end
    temp = (std(car_judge(1:SUM_BG))).^2;
    
     Z(p)=temp;
    if temp > threshold
        disp(['第',num2str(p),'个停车位：']);
    else
        disp(['第',num2str(p),'个停车位：']);
    end
    disp(['标定后方差为',num2str(temp)]);
    SUM_BG = 0;
    SUM_CAR = 0;
end

%% 占比确定

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
A_1 = (imread('2.jpg'));
A_1 = rgb2gray(A_1);
A_2=A_1 - A;
% A = rgb2gray(A);
% A_1_1 = im2bw(A,0.7);
% A_2_2 = im2bw(rgb2gray(A_2),0.7);       %有车图的二值化图
% %A_minus = A_2_2 - A_1_1;
% A_minus = abs(rgb2gray(A_2) - A);            %灰度值差
% A_minus2 = A_minus;
% %thresh = graythresh(A_minus);
% A_minus = im2bw(A_minus,0);


%% 显示差影运算
%figure;

%imshow(A_2);
%title('差影运算');
[m,n]=size(A_2);

w=fspecial('gaussian',[5 5],1.6);
A_2=imfilter(A_2,w,'replicate');

Q =[0,0,0];
z=0;
for i = (470+compensation_x):(590+compensation_y)
    for j = (12+compensation_x):(280)
        if A_2(i,j)>0
            z=z+1;
        end
    end
end
q=double(z/(m*n));
Q(1)=q;
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
Q(2)=q;
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
Q(3)=q;
 disp(['第三个车位的比例为',num2str(q)]);
 
 
 %% 判断条件/贝叶斯决策
 %% 数据读取

[n,t,r]=xlsread('样本1.xlsx');
l=1;m=1;
for k=1:540     %样本容量为540
  if n(k,3)==1
    you_1(l)=n(k,1);
    you_2(l)=n(k,2);
    l=l+1;
  else
    wu_1(m)=n(k,1);
    wu_2(m)=n(k,2);
    m=m+1;
  end
end

%% 车位有无车方差的直方图对比
%figure(1)
%subplot(2,1,1),
%hist(you_1), title('车位有车方差直方图'),xlabel('方差'),ylabel('频数');
%subplot(2,1,2),
%hist(wu_1), title('车位无车方差直方图'),xlabel('方差'),ylabel('频数');

%% 采用最大似然估计方法，求车位有无车方差比例分布的参数
% 假定车位有无车比例服从正态分布，此时极大似然估计对均值的估计结果为样本均值，即
you_1_u=mean(you_1);wu_1_u=mean(wu_1);        %方差均值
you_2_u=mean(you_2);you_2_u=mean(wu_2);        %比例均值 
you_1_S=0;you_2_S=0;wu_1_S=0;wu_2_S=0;
for m1=1:l-1
    you_1_S=you_1_S+(you_1(m1)-you_1_u)^2;   %有车方差与均值差值平方和
    you_2_S=you_2_S+(you_2(m1)-you_2_u)^2;   %有车比例与均值差值平方和
end
for w1=1:m-1
    wu_1_S=wu_1_S+(wu_1(w1)-wu_1_u)^2;   %无车方差与均值差值平方和
    wu_2_S=wu_2_S+(wu_2(w1)-you_2_u)^2;   %无车比例与均值差值平方和
end 
%极大似然估计对方差的估计结果为样本方差，即
sigma1=you_1_S/(l-1);
sigma2=you_2_S/(l-1);
sigma3=wu_1_S/(m-1);
sigma4=wu_2_S/(m-1);  
fprintf('有车方差极大似然估计均值为%f,方差为%f\n',you_1_u,sigma1);
fprintf('有车比例极大似然估计均值为%f,方差为%f\n',you_2_u,sigma2);
fprintf('无车方差极大似然估计均值为%f,方差为%f\n',wu_1_u,sigma3);
fprintf('无车比例极大似然估计均值为%f,方差为%f\n',you_2_u,sigma4); 

%% 先假定概率密度服从极大似然估计量所求分布，即 p(u1)~N(man_h_u,sigma1)...
%设置均值u的先验分布的方差均为10
sigmaX1=1206671;    
sigmaX2=0.00042;
sigmaX3=267.699233;
sigmaX4=0.000013;

%下面采用贝叶斯估计方法求有无车方差集体中分布的 u=uN

%由于p(u/X)~N(uN,sigmaN)
%最小错误率贝叶斯估计所估计出的均值为
uN1=(sigma1*sum(you_1)+sigmaX1*you_1_u)/(sigma1*length(you_1)+sigmaX1);
uN2=(sigma2*sum(you_2)+sigmaX2*you_2_u)/(sigma2*length(you_2)+sigmaX2);
uN3=(sigma3*sum(wu_1)+sigmaX3*wu_1_u)/(sigma3*length(wu_1)+sigmaX3);
uN4=(sigma4*sum(wu_2)+sigmaX4*you_2_u)/(sigma4*length(wu_2)+sigmaX4);

fprintf('\n\n\n');
fprintf('下面采用贝叶斯估计方法，求有无车方差和比例分布的参数\n');
fprintf('设置均值u的先验分布的方差均为10，即sigmaX1=10，sigmaX2=10，sigmaX3=10，sigmaX4=10\n');
fprintf('最小错误率贝叶斯估计所估计出的均值为%f,%f,%f,%f\n',uN1,uN2,uN3,uN4);
%% 采用最小错误率贝叶斯决策，画出类别判定的决策面
%协方差矩阵的计算过程,主对角元素即为sigma1,sigma2及sigma3,sigma4,
%只需求次对角线元素sigma12,sigma21及sigma34,sigma43，并且次对角线元素值相等sigma12=sigma21,sigma34=sigma43
%这里系数取N
C12=0;C34=0;      
for m2=1:length(you_1)
    C12=C12+(you_1(m2)-you_1_u)*(you_2(m2)-you_2_u);   
end
for w2=1:length(wu_1)
    C34=C34+(wu_1(w2)-wu_1_u)*(wu_2(w2)-you_2_u); 
end
sigma12=C12/m2;
sigma34=C34/w2;
%于是所求协方差矩阵为
sigma_you=[sigma1,sigma12;sigma12,sigma2];   
sigma_wu=[sigma3,sigma34;sigma34,sigma4]; 
%且易知先验概率为
p_you=m1/(m1+w1);   
p_wu=1-p_you;
syms  x1 x2 %定义矩阵x的两个元素
N_1=[x1-you_1_u,x2-you_2_u];    %(x-u1)
N_2=[x1-wu_1_u,x2-you_2_u];%(x-u2) 
%得决策面方程为
g=0.5*N_1*(sigma_you^(-1))*N_1'-0.5*N_2*(sigma_wu^(-1))*N_2'+0.5*log(det(sigma_you)/det(sigma_wu))-log(p_you/p_wu);
%图表2—最小错误率贝叶斯决策类别判定决策面图，及对样本判断。
figure(2)
%图表2：决策面
h=ezplot(g,[0,500,0,0.03]);title('决策面函数图像'),xlabel('方差'),ylabel('比例');
set(h,'Color','r') %将决策面边界标为红色
hold on;
%标出样本中有车方差比例数据（蓝色点）
for m3=1:l-1
    x=you_1(m3);y=you_2(m3);plot(x,y,'.');
end
%标出样本中无车方差比例数据（洋红色点）
for w3=1:m-1
    x=wu_1(w3);y=wu_2(w3);plot(x,y,'m.');
end




 %第一个停车位方差比例
 x=Z(1);y=Q(1);
N_1=[x-you_1_u,y-you_2_u];    %(x-u1)
N_2=[x-wu_1_u,y-you_2_u];%(x-u2) 

 plot(x,y,'c*');
if  0.5*N_1*(sigma_you^(-1))*N_1'-0.5*N_2*(sigma_wu^(-1))*N_2'+0.5*log(det(sigma_you)/det(sigma_wu))-log(p_you/p_wu)<0
        disp('第一个停车位：有车');
else
        disp('第一个停车位：无车');
end

 %第二个停车位方差比例
  x=Z(2);y=Q(2);
  
  N_1=[x-you_1_u,y-you_2_u];    %(x-u1)
N_2=[x-wu_1_u,y-you_2_u];%(x-u2) 
  plot(x,y,'g*'); 
if  0.5*N_1*(sigma_you^(-1))*N_1'-0.5*N_2*(sigma_wu^(-1))*N_2'+0.5*log(det(sigma_you)/det(sigma_wu))-log(p_you/p_wu)<0
    disp('第二个停车位：有车');
else
    disp('第二个停车位：无车');
end
  
  
  %第三个停车位方差比例
  x=Z(3);y=Q(3);
  N_1=[x-you_1_u,y-you_2_u];    %(x-u1)
N_2=[x-wu_1_u,y-you_2_u];%(x-u2) 
  plot(x,y,'b*'); 
if  0.5*N_1*(sigma_you^(-1))*N_1'-0.5*N_2*(sigma_wu^(-1))*N_2'+0.5*log(det(sigma_you)/det(sigma_wu))-log(p_you/p_wu)<0
     disp('第三个停车位：有车');
else
     disp('第三个停车位：无车');
end
  
  
  hold off;
  grid on;
 
 
