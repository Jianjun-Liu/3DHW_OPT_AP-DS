function [c,ceq] = st2(x)
      %%%%此程序用来计算套管位置c1、c2、c3(与文献SPE79164进行对比)  圆柱螺线法
 format long;
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
 k3 = x(11);  %第二圆弧段造斜率
 k5 = x(12);  %第三圆弧段造斜率
 Dkop = x(13);  %第一造斜点垂深
 Dd = x(14);  %第二造斜点垂深
 Db = x(15);  %第三造斜点垂深
 %以上为待优化的参数
 
%  a1 = 10;  %第一稳斜段井斜角
%  a2 = 40;  %第二稳斜段井斜角
%  a3 = 90;
%  fi1 = 270;
%  fi2 = 280;
%  fi3 = 270;
%  fi4 = 340;
%  fi5 = 340;
%  fi6 = 355;
%  k1 = 0.829;  %第一圆弧段造斜率
%  k3 = 0.91;  %第二圆弧段造斜率
%  k5 = 3.01;  %第三圆弧段造斜率
%  Dkop = 1000;  %第一造斜点垂深
%  Dd = 7000;  %第二造斜点垂深
%  Db = 10199.92;  %第三造斜点垂深
 
 
  a1 = a1 *pi/180;
  a2 = a2 *pi/180;
  a3 = a3 *pi/180;
  fi1 = fi1 * pi/180;
  fi2 = fi2 * pi/180;
  fi3 = fi3 * pi/180;
  fi4 = fi4 * pi/180;
  fi5 = fi5 * pi/180;
  fi6 = fi6 * pi/180;

R1 = 18000/(pi * k1);
R3 = 18000/(pi * k3);
R5 = 18000/(pi * k5);
D1 = R1 * ((fi2 - fi1)^2 * (sin(a1/2))^4 + (a1)^2)^0.5;%第一圆弧段长度
D3 = R3 * ((fi4 - fi3)^2 * (sin((a1 + a2)/2))^4 + (a2 - a1)^2)^0.5;
D5 = R5 * ((fi6 - fi5)^2 * (sin((a2 + a3)/2))^4 + (a3 - a2)^2)^0.5;

%D1 = R1 * ((fi2 - fi1)^2 * (sin(a1/2))^2 + (a1)^2)^0.5;
%D3 = R3 * ((fi4 - fi3)^2 * (sin((a1 + a2)/2))^2 + (a2 - a1)^2)^0.5;
%D5 = R5 * ((fi6 - fi5)^2 * (sin((a2 + a3)/2))^2 + (a3 - a2)^2)^0.5;

D2 = (Dd - Dkop - D1 * sin(a1)/a1) / cos(a1);%第一稳斜段长度
D4 = (Db - Dd - D3 * (sin(a2) - sin(a1))/(a2 - a1))/cos(a2);%第二稳斜段长度
c1 = Dkop + D1 *sin(a1)/a1;
%c2 = Dd + D3 * (sin(a2) - sin(a1))/(a2 - a1);
c2 = Dkop +  D1 * sin(a1)/a1 + D2 * cos(a1) + D3 * (sin(a2) - sin(a1))/(a2 - a1);
c3 = Dkop +  D1 * sin(a1)/a1 + D2 * cos(a1) + D3 * (sin(a2) - sin(a1))/(a2 - a1) + D4 *cos(a2) +D5 * (sin(a3) - sin(a2))/(a3 - a2);
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
c=[ -D2 ; -D4 ; (fi3 - fi2) * 180/pi/D2 - 5/100; (fi5 - fi4) * 180/pi/D4 - 5/100; 1800 - c1; c1 - 2200; 7200 - c2;  c2 - 8700 ; 10850 - c3; c3 - 10900  ];
ceq=[];

end
