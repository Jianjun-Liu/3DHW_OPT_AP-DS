function [c,ceq] = st1(x)
      %%%%�˳������������׹�λ��c1��c2��c3(������SPE79164���жԱ�)  б��Բ����
 format long;
 %��һбֱ�߶ξ�б��
 a1 = x(1);  %��һ��б�ξ�б��
 a2 = x(2);  %�ڶ���б�ξ�б��
 a3 = x(3);
 fi1 = x(4);
 fi2 = x(5);
 fi3 = x(6);
 fi4 = x(7);
 fi5 = x(8);
 fi6 = x(9);
 k1 = x(10);  %��һԲ������б��
 k2 = x(11);  %�ڶ�Բ������б��
 k3 = x(12);  %����Բ������б��
 Dkop = x(13);  %��һ��б�㴹��
 Dd = x(14);  %�ڶ���б�㴹��
 Db = x(15);  %������б�㴹��
 %����Ϊ���Ż��Ĳ���
 a0 = 0;
% a1=9.95;
% fi1=50;
% a2=40.73;
% fi2=67.40;
% k1=6;
% k2=6.786;
% Dkop=1121.31;
% Dc=1850.8;

r1 = acos(cosd(a1) * cosd(a0) + sind(a1) * sind(a0) * cosd(fi2 - fi1));%��һ��б�Σ�Բ���Σ����Ƚ�
%r1 = 2 * asin(sqrt(sind( (a1- a0)/2 )^2 + sind( (fi2- fi1)/2 )^2 * sind(a0) * sind(a1) ));

r2 = acos(cosd(a2) * cosd(a1) + sind(a2) * sind(a1) * cosd(fi4 - fi3));%�ڶ���б�Σ�Բ���Σ����Ƚ�
%r2 = 2 * asin(sqrt(sind( (a2- a1)/2 )^2 + sind( (fi4- fi3)/2 )^2 * sind(a1) * sind(a2) ));

r3 = acos(cosd(a2) * cosd(a3) + sind(a2) * sind(a3) * cosd(fi6 - fi5));%�ڶ���б�Σ�Բ���Σ����Ƚ�
%r3 = 2 * asin(sqrt(sind( (a3- a2)/2 )^2 + sind( (fi6- fi5)/2 )^2 * sind(a2) * sind(a3) )); 

D2 = (Dd - Dkop - 100/(k1 * pi/180) * tan(r1/2) * (cosd(a1) + cosd(a0)))/cosd(a1);%��һ��б�γ���
D4 = (Db - Dd - 100/(k2 * pi/180) * tan(r2/2) * (cosd(a2) + cosd(a1)))/cosd(a2);%��һ��б�γ���

c1 = Dkop +  100/(k1 * pi/180) * tan(r1/2) * (cosd(a1) + cosd(a0));
c2 = Dd + 100/(k2 * pi/180) * tan(r2/2) * (cosd(a2) + cosd(a1));
c3 = Db + 100/(k3 * pi/180) * tan(r3/2) * (cosd(a3) + cosd(a2));
%��c3 = TVD��

% %����N,E,H�ֱ��ʾбֱ�߶κ�Բ���ε�N,E,H���������ļ�����ʽ
% N1=Dkop*tand(a1)*cosd(fi1);
% E1=Dkop*tand(a1)*sind(fi1);
% N2=30/(k1*pi/180)*tan(r1/2)*(sind(a1)*cosd(fi1)+sind(a2)*cosd(fi2));
% E2=30/(k1*pi/180)*tan(r1/2)*(sind(a1)*sind(fi1)+sind(a2)*sind(fi2));
% N3=D2*sind(a2)*cosd(fi2);
% E3=D2*sind(a2)*sind(fi2);
% 
% N4=30/(k2*pi/180)*tan(r2/2)*(sind(a2)*cosd(fi2)+sind(a3)*cosd(fi3));
% E4=30/(k2*pi/180)*tan(r2/2)*(sind(a2)*sind(fi2)+sind(a3)*sind(fi3));
% H4=30/(k2*pi/180)*tan(r2/2)*(cosd(a2)+cosd(a3));
% %����Ϊ�ܵ�N,E,H�����������ʽ����ֵӦ���ڸ������װе�����
% N=N1+N2+N3+N4;
% E=E1+E2+E3+E4;
% H=H4+Dc;
c=[ -D2 ; -D4 ; (fi3 - fi2)/D2 - 5/100; (fi5 - fi4)/D4 - 5/100; 1800 - c1; c1 - 2200;  7200 - c2;  c2 - 8700 ; 10850 - c3; c3 - 10900];
ceq=[];

end
