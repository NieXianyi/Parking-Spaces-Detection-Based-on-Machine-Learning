
clc;
close all;
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
figure(1)
subplot(2,1,1),
hist(you_1), title('��λ�г�����ֱ��ͼ'),xlabel('����'),ylabel('Ƶ��');
subplot(2,1,2),
hist(wu_1), title('��λ�޳�����ֱ��ͼ'),xlabel('����'),ylabel('Ƶ��');

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
% %�ж������ķ�������ֱ�Ϊ(160,45)
 x=3000;y=0.019;
 plot(x,y,'c*');
 %�ж������ķ�������ֱ�Ϊ(7.012,0.002)
  x=7.012;y=0.002;
  plot(x,y,'g*'); 
  hold off;
  grid on;

