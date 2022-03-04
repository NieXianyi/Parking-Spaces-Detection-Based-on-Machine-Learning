clear;
close all;
%% �û�����
H = 2; %������߶�
L = 4; %��������
l = 2;  %���Ŀ��
h = 1.5; %���ĸ߶�


up=[195,435;
    338,435;
    481,435;
   618,435];
down=[1,542;
    278,540;
    548,538;
    805,536];



threshold = 20; %���������ռ��λ����İٷֱȵ���ֵ
%threshold2 = 15; %���������ռ��λ����İٷֱȵ���ֵ
%% ��ȡͼ���Լ�Ԥ����
%A = (imread('�ճ�λ(��չ)(���ȵ���).bmp'));
A = (imread('����.jpg'));
A_2 = (imread('6.jpg'));
A = rgb2gray(A);
A_1_1 = im2bw(A,0.7);
A_2_2 = im2bw(rgb2gray(A_2),0.7);       %�г�ͼ�Ķ�ֵ��ͼ
%A_minus = A_2_2 - A_1_1;
A_minus = abs(rgb2gray(A_2) - A);            %�Ҷ�ֵ��
A_minus2 = A_minus;
%thresh = graythresh(A_minus);
A_minus = im2bw(A_minus,0);

subplot(221);
imshow(A);
title('ͣ����');
subplot(222);
imshow(A_2);


title('��ͣ�ų�����ͣ����');
subplot(2,2,4);
imshow(A_minus2);
title('��Ӱ����');
%% ������궨
theta1 = atan((H-h)/L);
demarcate_l = sin(theta1)*h/l; 
%�õ��ӽ�����������ϵ�ı���

compensation_x = 0;         %һ�㲻��Ҫ�궨x����
compensation_y = -round(demarcate_l*150);  %��Ҫ��Y���򲹳����ɸ�����
%% �����������λ����
for p = 1:1:size(down,1)
    for i = 1:1:down(p,2)-up(p,2)+1
        k = (up(p,1)-down(p,1))/(up(p,2)-down(p,2));
        B{p}(i,1) = round(up(p,1)+i*k);
        B{p}(i,2) = up(p,2)+i-1;
    end
end
    
%% ��ԭͼ��ǳ���λ
 subplot(223);
 imshow(A);
 title('���ͣ��λ');
 hold on
 for k = 1:length(B)
     %plot(B{k}(:,1), B{k}(:,2), 'G', 'LineWidth', 2)
 end
 for k = 1:length(B)-1
     %plot([B{k}(1,1),B{k+1}(1,1)],[B{k}(1,2),B{k+1}(1,2)], 'G', 'LineWidth', 2);
    % plot([B{k}(size(B{k},1),1),B{k+1}(size(B{k+1},1),1)],[B{k}(size(B{k},1),2),B{k+1}(size(B{k+1},1),2)], 'G', 'LineWidth', 2);
 end
%% ��ʾ��Ӱ����
figure;
subplot(221);
imshow(A_minus);
title('��Ӱ����');
%% ������
A_minus = bwmorph(A_minus,'close');% ������
%% ��ʾ������
subplot(2,2,2);
imshow(A_minus);
title('������');
% %% �ڱ�����ͼ���ͣ��λ
subplot(2,2,3);
imshow(A_minus);
title('ͣ��λ���');
hold on
for k = 1:length(B)
   % plot(B{k}(:,1), B{k}(:,2), 'G', 'LineWidth', 2)
end
for k = 1:length(B)-1
   % plot([B{k}(1,1),B{k+1}(1,1)],[B{k}(1,2),B{k+1}(1,2)], 'G', 'LineWidth', 2);
  %  plot([B{k}(size(B{k},1),1),B{k+1}(size(B{k+1},1),1)],[B{k}(size(B{k},1),2),B{k+1}(size(B{k+1},1),2)], 'G', 'LineWidth', 2);
end
%% ɨ��ÿ�����ؼ���Ƿ��г�
subplot(2,2,4);
imshow(A_2);
title('��ǽ��');%������Ϊ�ڵ�������δ�ڵ����ִ�С��ͬ
hold on
SUM_BG = 0;%SUM_CAR = 0;%�����복���ص�
 car_judge = [];
for p = 1:1:length(B)-1          %����û�г�
    B1 = max(B{p}(:,2)) - min(B{p}(:,2));
    B2 = max(B{p+1}(:,2)) - min(B{p+1}(:,2));   %Ѱ�������Сֵ֮��
    for i = 1:1:min(B1,B2)        %��һѭ������С��ֵ
        for j = B{p}(i,1):1:B{p+1}(i,1)     %����߱������ұ�
            SUM_BG = SUM_BG + 1;            %���ظ���
            car_judge(SUM_BG) = A_minus2(B{p}(1,2)+i-1,j);
           % plot([j,j+1],[B{p}(1,2)+i-1,B{p}(1,2)+i-1], 'B', 'LineWidth', 2);
            if A_minus(B{p}(1,2)+i-1,j)== 1   %ʶ��
%                 SUM_CAR = SUM_CAR + 1;
                %plot([j,j+1],[B{p}(1,2)+i-1,B{p}(1,2)+i-1], 'R', 'LineWidth', 2);                
            end
        end
    end
    temp = (std(car_judge(1:SUM_BG))).^2;
    if temp > threshold
        disp(['δ�궨ʱ��',num2str(p),'��ͣ��λ���г�']);
    else
        disp(['δ�궨ʱ��',num2str(p),'��ͣ��λ���޳�']);
    end
    disp(['δ�궨ʱ����Ϊ',num2str(temp)]);
    SUM_BG = 0;
end

%% �ָ���

%������������������������������������������������������������������������������������%

%% ��ʾ�궨��ͣ��λʶ������
figure
 subplot(221);
 imshow(A);
title('�궨��ͣ��λ���(�������)');
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
% ��ʾ�궨���ͣ��λ����
subplot(222);
imshow(A);
title('�궨���ͣ��λ���(�������)');
hold on
for k = 1:length(B)
    %plot(B{k}(:,1)+compensation_x, B{k}(:,2)+compensation_y, 'G', 'LineWidth', 2)
end
for k = 1:length(B)-1
    %plot([B{k}(1,1)+compensation_x,B{k+1}(1,1)+compensation_x],[B{k}(1,2)+compensation_y,B{k+1}(1,2)+compensation_y], 'G',[1,0,1], 2);
   % plot([B{k}(size(B{k},1)+compensation_x,1),B{k+1}(size(B{k+1},1)+compensation_x,1)],[B{k}(size(B{k},1)+compensation_y,2),B{k+1}(size(B{k+1},1)+compensation_y,2)], 'G',[1,0,1], 2);
end
% ��ʾ�궨�����ʾ���
subplot(2,2,3.5);
imshow(A_2);
title('�궨���ǽ��');%������Ϊ�ڵ�������δ�ڵ����ִ�С��ͬ
hold on
SUM_BG = 0;SUM_CAR = 0;%�����복���ص�

for p = 1:1:length(B)-1          %����û�г�
    B1 = max(B{p}(:,2)) - min(B{p}(:,2));
    B2 = max(B{p+1}(:,2)) - min(B{p+1}(:,2));   %Ѱ�������Сֵ֮��
    for i = 1:1:min(B1,B2)        %��һѭ������С��ֵ
        for j = B{p}(i,1):1:B{p+1}(i,1)     %����߱������ұ�
            SUM_BG = SUM_BG + 1;
            car_judge(SUM_BG) = A_minus2(B{p}(1,2)+i-1+compensation_y,j+compensation_x);
            %plot([j+compensation_x,j+1+compensation_x],[B{p}(1,2)+i-1+compensation_y,B{p}(1,2)+i-1+compensation_y], 'B', 'LineWidth', 2);
            if A_minus(B{p}(1,2)+i-1+compensation_y,j+compensation_x) == 1   %���벹��
                SUM_CAR = SUM_CAR + 1;
              %  plot([j+compensation_x,j+1+compensation_x],[B{p}(1,2)+i-1+compensation_y,B{p}(1,2)+i-1+compensation_y], 'R', 'LineWidth', 2);                
            end
        end
    end
    temp = (std(car_judge(1:SUM_BG))).^2;
    if temp > threshold
        disp(['��',num2str(p),'��ͣ��λ���г�']);
    else
        disp(['��',num2str(p),'��ͣ��λ���޳�']);
    end
    disp(['�궨�󷽲�Ϊ',num2str(temp)]);
    SUM_BG = 0;
    SUM_CAR = 0;
end