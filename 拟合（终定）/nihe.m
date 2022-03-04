
clc;
close all;
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
figure(1)
subplot(2,1,1),
hist(you_1), title('车位有车方差直方图'),xlabel('方差'),ylabel('频数');
subplot(2,1,2),
hist(wu_1), title('车位无车方差直方图'),xlabel('方差'),ylabel('频数');

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
% %判断样本的方差比例分别为(160,45)
 x=3000;y=0.019;
 plot(x,y,'c*');
 %判断样本的方差比例分别为(7.012,0.002)
  x=7.012;y=0.002;
  plot(x,y,'g*'); 
  hold off;
  grid on;

