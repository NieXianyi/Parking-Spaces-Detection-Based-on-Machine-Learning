%% ��һ��

clc;
close all;
%% �û�����
H = 15; %������߶�
L = 30; %��������
l = 2;  %���Ŀ��
h = 1.5; %���ĸ߶�

up = [210,470;          %�������ɵ�
    352,469;
    502,468;
    652,466];
down = [12,588;        %�������ɵ�
    280,589;
    560,589;
    830,590];

%threshold = 20; %���������ռ��λ����İٷֱȵ���ֵ
%threshold2 = 15; %���������ռ��λ����İٷֱȵ���ֵ

%% ������궨
theta1 = atan((H-h)/L);
demarcate_l = sin(theta1)*h/l; 
%�õ��ӽ�����������ϵ�ı���

compensation_x = 0;         %һ�㲻��Ҫ�궨x����
compensation_y = round(demarcate_l*150);  %��Ҫ��Y���򲹳����ɸ�����

%% ��ȡͼ���Լ�Ԥ����
%A = (imread('�ճ�λ(��չ)(���ȵ���).bmp'));
A = (imread('û�г�.jpg'));
A = rgb2gray(A);
A_1 = (imread('12.jpg'));
A_1 = rgb2gray(A_1);
A_2=A_1-A;
% A = rgb2gray(A);
% A_1_1 = im2bw(A,0.7);
% A_2_2 = im2bw(rgb2gray(A_2),0.7);       %�г�ͼ�Ķ�ֵ��ͼ
% %A_minus = A_2_2 - A_1_1;
% A_minus = abs(rgb2gray(A_2) - A);            %�Ҷ�ֵ��
% A_minus2 = A_minus;
% %thresh = graythresh(A_minus);
% A_minus = im2bw(A_minus,0);


%% ��ʾ��Ӱ����
figure;

imshow(A_2);
title('��Ӱ����');
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

 disp(['��һ����λ�ı���Ϊ',num2str(q)]);
 
 z=0;
for i = (469+compensation_x):(589+compensation_y)
    for j = (352+compensation_x):(502)
        if A_2(i,j)>0
            z=z+1;
        end
    end
end
q=double(z/(m*n));
 disp(['�ڶ�����λ�ı���Ϊ',num2str(q)]);
 
  z=0;
for i = (468+compensation_x):(589+compensation_y)
    for j = (560+compensation_x):(830)
        if A_2(i,j)>0
            z=z+1;
        end
    end
end
q=double(z/(m*n));
 disp(['��������λ�ı���Ϊ',num2str(q)]);