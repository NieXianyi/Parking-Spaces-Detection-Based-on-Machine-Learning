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
%% �û�����
H = 2; %������߶�
L = 4; %��������
l = 2;  %���Ŀ��
h = 1.5; %���ĸ߶�

up = [205,467;          %�������ɵ�
    350,470;
    505,470;
    655,470];
down = [10,590;        %�������ɵ�
    280,590;
    556,590;
    850,590];

threshold = 20; %���������ռ��λ����İٷֱȵ���ֵ
%threshold2 = 15; %���������ռ��λ����İٷֱȵ���ֵ
%% ��ȡͼ���Լ�Ԥ����
%A = (imread('�ճ�λ(��չ)(���ȵ���).bmp'));
A = (imread('û�г�.jpg'));
A_2 = (imread('2.jpg'));
A = rgb2gray(A);
A_1_1 = im2bw(A,0.7);
A_2_2 = im2bw(rgb2gray(A_2),0.7);       %�г�ͼ�Ķ�ֵ��ͼ
%A_minus = A_2_2 - A_1_1;
A_minus = abs(rgb2gray(A_2) - A);            %�Ҷ�ֵ��
A_minus2 = A_minus;
%thresh = graythresh(A_minus);
A_minus = im2bw(A_minus,0);

% subplot(221);
% imshow(A);
% title('ͣ����');
% subplot(222);
% imshow(A_2);


% title('��ͣ�ų�����ͣ����');
% subplot(2,2,4);
% imshow(A_minus2);
% title('��Ӱ����');
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
%  subplot(223);
%  imshow(A);
%  title('���ͣ��λ');
 hold on
 for k = 1:length(B)
    plot(B{k}(:,1), B{k}(:,2), 'G', 'LineWidth', 2)
 end
 for k = 1:length(B)-1
     %plot([B{k}(1,1),B{k+1}(1,1)],[B{k}(1,2),B{k+1}(1,2)], 'G', 'LineWidth', 2);
     %plot([B{k}(size(B{k},1),1),B{k+1}(size(B{k+1},1),1)],[B{k}(size(B{k},1),2),B{k+1}(size(B{k+1},1),2)], 'G', 'LineWidth', 2);
 end
%% ��ʾ��Ӱ����
% figure;
% subplot(221);
% imshow(A_minus);
% title('��Ӱ����');
%% ������
A_minus = bwmorph(A_minus,'close');% ������
%% ��ʾ������
 %subplot(2,2,2);
 %imshow(A_minus);
 %title('������');
 %% �ڱ�����ͼ���ͣ��λ
 %subplot(2,2,3);
 %imshow(A_minus);
 %title('ͣ��λ���');
hold on
for k = 1:length(B)
   % plot(B{k}(:,1), B{k}(:,2), 'G', 'LineWidth', 2)
end
for k = 1:length(B)-1
    %plot([B{k}(1,1),B{k+1}(1,1)],[B{k}(1,2),B{k+1}(1,2)], 'G', 'LineWidth', 2);
    %plot([B{k}(size(B{k},1),1),B{k+1}(size(B{k+1},1),1)],[B{k}(size(B{k},1),2),B{k+1}(size(B{k+1},1),2)], 'G', 'LineWidth', 2);
end
%% ɨ��ÿ�����ؼ���Ƿ��г�
% subplot(2,2,4);
% imshow(A_2);
% title('��ǽ��');%������Ϊ�ڵ�������δ�ڵ����ִ�С��ͬ
hold on
SUM_BG = 0;SUM_CAR = 0;%�����복���ص�
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
                 SUM_CAR = SUM_CAR + 1;
               % plot([j,j+1],[B{p}(1,2)+i-1,B{p}(1,2)+i-1], 'R', 'LineWidth', 2);                
            end
        end
    end
    temp = (std(car_judge(1:SUM_BG))).^2;
    SUM_BG = 0;
end

%% ��ʾ�궨��ͣ��λʶ������
% figure
%  subplot(221);
%  imshow(A);
% title('�궨��ͣ��λ���(�������)');
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
% ��ʾ�궨���ͣ��λ����
% subplot(222);
% imshow(A);
% title('�궨���ͣ��λ���(�������)');
hold on
for k = 1:length(B)
    plot(B{k}(:,1)+compensation_x, B{k}(:,2)+compensation_y, 'G', 'LineWidth', 2)
end
for k = 1:length(B)-1
    plot([B{k}(1,1)+compensation_x,B{k+1}(1,1)+compensation_x],[B{k}(1,2)+compensation_y,B{k+1}(1,2)+compensation_y], 'G',[1,0,1], 2);
    plot([B{k}(size(B{k},1)+compensation_x,1),B{k+1}(size(B{k+1},1)+compensation_x,1)],[B{k}(size(B{k},1)+compensation_y,2),B{k+1}(size(B{k+1},1)+compensation_y,2)], 'G',[1,0,1], 2);
end
% ��ʾ�궨�����ʾ���
% subplot(2,2,3.5);
% imshow(A_2);
% title('�궨���ǽ��');%������Ϊ�ڵ�������δ�ڵ����ִ�С��ͬ
hold on

SUM_BG = 0;SUM_CAR = 0;%�����복���ص�
%����û�г�

Z = [0,0,0];
for p = 1:1:length(B)-1          
    B1 = max(B{p}(:,2)) - min(B{p}(:,2));
    B2 = max(B{p+1}(:,2)) - min(B{p+1}(:,2));   %Ѱ�������Сֵ֮��
    for i = 1:1:min(B1,B2)        %��һѭ������С��ֵ
        for j = B{p}(i,1):1:B{p+1}(i,1)     %����߱������ұ�
            SUM_BG = SUM_BG + 1;
            car_judge(SUM_BG) = A_minus2(B{p}(1,2)+i-1+compensation_y,j+compensation_x);
          % plot([j+compensation_x,j+1+compensation_x],[B{p}(1,2)+i-1+compensation_y,B{p}(1,2)+i-1+compensation_y], 'B', 'LineWidth', 2);
           
          if A_minus(B{p}(1,2)+i-1+compensation_y,j+compensation_x) == 1   %���벹��
                SUM_CAR = SUM_CAR + 1;
             %plot([j+compensation_x,j+1+compensation_x],[B{p}(1,2)+i-1+compensation_y,B{p}(1,2)+i-1+compensation_y], 'R', 'LineWidth', 2);                
          end
            
        end
    end
    temp = (std(car_judge(1:SUM_BG))).^2;
    
     Z(p)=temp;
    if temp > threshold
        disp(['��',num2str(p),'��ͣ��λ��']);
    else
        disp(['��',num2str(p),'��ͣ��λ��']);
    end
    disp(['�궨�󷽲�Ϊ',num2str(temp)]);
    SUM_BG = 0;
    SUM_CAR = 0;
end

%% ռ��ȷ��

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
A_1 = (imread('2.jpg'));
A_1 = rgb2gray(A_1);
A_2=A_1 - A;
% A = rgb2gray(A);
% A_1_1 = im2bw(A,0.7);
% A_2_2 = im2bw(rgb2gray(A_2),0.7);       %�г�ͼ�Ķ�ֵ��ͼ
% %A_minus = A_2_2 - A_1_1;
% A_minus = abs(rgb2gray(A_2) - A);            %�Ҷ�ֵ��
% A_minus2 = A_minus;
% %thresh = graythresh(A_minus);
% A_minus = im2bw(A_minus,0);


%% ��ʾ��Ӱ����
%figure;

%imshow(A_2);
%title('��Ӱ����');
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
Q(2)=q;
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
Q(3)=q;
 disp(['��������λ�ı���Ϊ',num2str(q)]);
 
 
 %% �ж�����/��Ҷ˹����
 %% ���ݶ�ȡ

[n,t,r]=xlsread('����1.xlsx');
l=1;m=1;
for k=1:540     %��������Ϊ540
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

%% ��λ���޳������ֱ��ͼ�Ա�
%figure(1)
%subplot(2,1,1),
%hist(you_1), title('��λ�г�����ֱ��ͼ'),xlabel('����'),ylabel('Ƶ��');
%subplot(2,1,2),
%hist(wu_1), title('��λ�޳�����ֱ��ͼ'),xlabel('����'),ylabel('Ƶ��');

%% ���������Ȼ���Ʒ�������λ���޳���������ֲ��Ĳ���
% �ٶ���λ���޳�����������̬�ֲ�����ʱ������Ȼ���ƶԾ�ֵ�Ĺ��ƽ��Ϊ������ֵ����
you_1_u=mean(you_1);wu_1_u=mean(wu_1);        %�����ֵ
you_2_u=mean(you_2);you_2_u=mean(wu_2);        %������ֵ 
you_1_S=0;you_2_S=0;wu_1_S=0;wu_2_S=0;
for m1=1:l-1
    you_1_S=you_1_S+(you_1(m1)-you_1_u)^2;   %�г��������ֵ��ֵƽ����
    you_2_S=you_2_S+(you_2(m1)-you_2_u)^2;   %�г��������ֵ��ֵƽ����
end
for w1=1:m-1
    wu_1_S=wu_1_S+(wu_1(w1)-wu_1_u)^2;   %�޳��������ֵ��ֵƽ����
    wu_2_S=wu_2_S+(wu_2(w1)-you_2_u)^2;   %�޳��������ֵ��ֵƽ����
end 
%������Ȼ���ƶԷ���Ĺ��ƽ��Ϊ���������
sigma1=you_1_S/(l-1);
sigma2=you_2_S/(l-1);
sigma3=wu_1_S/(m-1);
sigma4=wu_2_S/(m-1);  
fprintf('�г��������Ȼ���ƾ�ֵΪ%f,����Ϊ%f\n',you_1_u,sigma1);
fprintf('�г�����������Ȼ���ƾ�ֵΪ%f,����Ϊ%f\n',you_2_u,sigma2);
fprintf('�޳��������Ȼ���ƾ�ֵΪ%f,����Ϊ%f\n',wu_1_u,sigma3);
fprintf('�޳�����������Ȼ���ƾ�ֵΪ%f,����Ϊ%f\n',you_2_u,sigma4); 

%% �ȼٶ������ܶȷ��Ӽ�����Ȼ����������ֲ����� p(u1)~N(man_h_u,sigma1)...
%���þ�ֵu������ֲ��ķ����Ϊ10
sigmaX1=1206671;    
sigmaX2=0.00042;
sigmaX3=267.699233;
sigmaX4=0.000013;

%������ñ�Ҷ˹���Ʒ��������޳�������зֲ��� u=uN

%����p(u/X)~N(uN,sigmaN)
%��С�����ʱ�Ҷ˹���������Ƴ��ľ�ֵΪ
uN1=(sigma1*sum(you_1)+sigmaX1*you_1_u)/(sigma1*length(you_1)+sigmaX1);
uN2=(sigma2*sum(you_2)+sigmaX2*you_2_u)/(sigma2*length(you_2)+sigmaX2);
uN3=(sigma3*sum(wu_1)+sigmaX3*wu_1_u)/(sigma3*length(wu_1)+sigmaX3);
uN4=(sigma4*sum(wu_2)+sigmaX4*you_2_u)/(sigma4*length(wu_2)+sigmaX4);

fprintf('\n\n\n');
fprintf('������ñ�Ҷ˹���Ʒ����������޳�����ͱ����ֲ��Ĳ���\n');
fprintf('���þ�ֵu������ֲ��ķ����Ϊ10����sigmaX1=10��sigmaX2=10��sigmaX3=10��sigmaX4=10\n');
fprintf('��С�����ʱ�Ҷ˹���������Ƴ��ľ�ֵΪ%f,%f,%f,%f\n',uN1,uN2,uN3,uN4);
%% ������С�����ʱ�Ҷ˹���ߣ���������ж��ľ�����
%Э�������ļ������,���Խ�Ԫ�ؼ�Ϊsigma1,sigma2��sigma3,sigma4,
%ֻ����ζԽ���Ԫ��sigma12,sigma21��sigma34,sigma43�����ҴζԽ���Ԫ��ֵ���sigma12=sigma21,sigma34=sigma43
%����ϵ��ȡN
C12=0;C34=0;      
for m2=1:length(you_1)
    C12=C12+(you_1(m2)-you_1_u)*(you_2(m2)-you_2_u);   
end
for w2=1:length(wu_1)
    C34=C34+(wu_1(w2)-wu_1_u)*(wu_2(w2)-you_2_u); 
end
sigma12=C12/m2;
sigma34=C34/w2;
%��������Э�������Ϊ
sigma_you=[sigma1,sigma12;sigma12,sigma2];   
sigma_wu=[sigma3,sigma34;sigma34,sigma4]; 
%����֪�������Ϊ
p_you=m1/(m1+w1);   
p_wu=1-p_you;
syms  x1 x2 %�������x������Ԫ��
N_1=[x1-you_1_u,x2-you_2_u];    %(x-u1)
N_2=[x1-wu_1_u,x2-you_2_u];%(x-u2) 
%�þ����淽��Ϊ
g=0.5*N_1*(sigma_you^(-1))*N_1'-0.5*N_2*(sigma_wu^(-1))*N_2'+0.5*log(det(sigma_you)/det(sigma_wu))-log(p_you/p_wu);
%ͼ��2����С�����ʱ�Ҷ˹��������ж�������ͼ�����������жϡ�
figure(2)
%ͼ��2��������
h=ezplot(g,[0,500,0,0.03]);title('�����溯��ͼ��'),xlabel('����'),ylabel('����');
set(h,'Color','r') %��������߽��Ϊ��ɫ
hold on;
%����������г�����������ݣ���ɫ�㣩
for m3=1:l-1
    x=you_1(m3);y=you_2(m3);plot(x,y,'.');
end
%����������޳�����������ݣ����ɫ�㣩
for w3=1:m-1
    x=wu_1(w3);y=wu_2(w3);plot(x,y,'m.');
end




 %��һ��ͣ��λ�������
 x=Z(1);y=Q(1);
N_1=[x-you_1_u,y-you_2_u];    %(x-u1)
N_2=[x-wu_1_u,y-you_2_u];%(x-u2) 

 plot(x,y,'c*');
if  0.5*N_1*(sigma_you^(-1))*N_1'-0.5*N_2*(sigma_wu^(-1))*N_2'+0.5*log(det(sigma_you)/det(sigma_wu))-log(p_you/p_wu)<0
        disp('��һ��ͣ��λ���г�');
else
        disp('��һ��ͣ��λ���޳�');
end

 %�ڶ���ͣ��λ�������
  x=Z(2);y=Q(2);
  
  N_1=[x-you_1_u,y-you_2_u];    %(x-u1)
N_2=[x-wu_1_u,y-you_2_u];%(x-u2) 
  plot(x,y,'g*'); 
if  0.5*N_1*(sigma_you^(-1))*N_1'-0.5*N_2*(sigma_wu^(-1))*N_2'+0.5*log(det(sigma_you)/det(sigma_wu))-log(p_you/p_wu)<0
    disp('�ڶ���ͣ��λ���г�');
else
    disp('�ڶ���ͣ��λ���޳�');
end
  
  
  %������ͣ��λ�������
  x=Z(3);y=Q(3);
  N_1=[x-you_1_u,y-you_2_u];    %(x-u1)
N_2=[x-wu_1_u,y-you_2_u];%(x-u2) 
  plot(x,y,'b*'); 
if  0.5*N_1*(sigma_you^(-1))*N_1'-0.5*N_2*(sigma_wu^(-1))*N_2'+0.5*log(det(sigma_you)/det(sigma_wu))-log(p_you/p_wu)<0
     disp('������ͣ��λ���г�');
else
     disp('������ͣ��λ���޳�');
end
  
  
  hold off;
  grid on;
 
 
