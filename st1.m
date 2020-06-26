function [c,ceq] = st1(x)
      %%%%此程序用来计算套管位置c1、c2、c3(与文献SPE79164进行对比)  斜面圆弧法
 format long;
 %第一斜直线段井斜角
 a1 = x(1);  %第一稳斜段井斜角
 a2 = x(2);  %第二稳斜段井斜角
 a3 = x(3);
 fi1 = x(4);
 fi2 = x(5);
 fi3 = x(6);
 fi4 = x(7);
 fi5 = x(8);
 fi6 = x(9);
 k1 = x(10);  %第一圆弧段造斜率
 k2 = x(11);  %第二圆弧段造斜率
 k3 = x(12);  %第三圆弧段造斜率
 Dkop = x(13);  %第一造斜点垂深
 Dd = x(14);  %第二造斜点垂深
 Db = x(15);  %第三造斜点垂深
 %以上为待优化的参数
 a0 = 0;
% a1=9.95;
% fi1=50;
% a2=40.73;
% fi2=67.40;
% k1=6;
% k2=6.786;
% Dkop=1121.31;
% Dc=1850.8;

r1 = acos(cosd(a1) * cosd(a0) + sind(a1) * sind(a0) * cosd(fi2 - fi1));%第一造斜段（圆弧段）狗腿角
%r1 = 2 * asin(sqrt(sind( (a1- a0)/2 )^2 + sind( (fi2- fi1)/2 )^2 * sind(a0) * sind(a1) ));

r2 = acos(cosd(a2) * cosd(a1) + sind(a2) * sind(a1) * cosd(fi4 - fi3));%第二造斜段（圆弧段）狗腿角
%r2 = 2 * asin(sqrt(sind( (a2- a1)/2 )^2 + sind( (fi4- fi3)/2 )^2 * sind(a1) * sind(a2) ));

r3 = acos(cosd(a2) * cosd(a3) + sind(a2) * sind(a3) * cosd(fi6 - fi5));%第二造斜段（圆弧段）狗腿角
%r3 = 2 * asin(sqrt(sind( (a3- a2)/2 )^2 + sind( (fi6- fi5)/2 )^2 * sind(a2) * sind(a3) )); 

D2 = (Dd - Dkop - 100/(k1 * pi/180) * tan(r1/2) * (cosd(a1) + cosd(a0)))/cosd(a1);%第一稳斜段长度
D4 = (Db - Dd - 100/(k2 * pi/180) * tan(r2/2) * (cosd(a2) + cosd(a1)))/cosd(a2);%第一稳斜段长度

c1 = Dkop +  100/(k1 * pi/180) * tan(r1/2) * (cosd(a1) + cosd(a0));
c2 = Dd + 100/(k2 * pi/180) * tan(r2/2) * (cosd(a2) + cosd(a1));
c3 = Db + 100/(k3 * pi/180) * tan(r3/2) * (cosd(a3) + cosd(a2));
%（c3 = TVD）

% %以下N,E,H分别表示斜直线段和圆弧段的N,E,H坐标增量的计算表达式
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
% %以下为总的N,E,H坐标增量表达式，其值应等于给定的首靶点坐标
% N=N1+N2+N3+N4;
% E=E1+E2+E3+E4;
% H=H4+Dc;
c=[ -D2 ; -D4 ; (fi3 - fi2)/D2 - 5/100; (fi5 - fi4)/D4 - 5/100; 1800 - c1; c1 - 2200;  7200 - c2;  c2 - 8700 ; 10850 - c3; c3 - 10900];
ceq=[];

end
