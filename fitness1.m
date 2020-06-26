function y=fitness1(x)
    %%%%此程序用来计算轨道总长度(与文献SPE79164进行对比)(迭代7000次L=14970)     斜面圆弧法
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
 k2 = x(11);  %第二圆弧段造斜率
 k3 = x(12);  %第三圆弧段造斜率
 Dkop = x(13);  %第一造斜点垂深
 Dd = x(14);  %第二造斜点垂深
 Db = x(15);  %第三造斜点垂深
 %以上为待优化的参数
 a0 = 0;

%  a1 = 10;  %第一稳斜段井斜角
%  a2 = 40;  %第二稳斜段井斜角
%  a3 = 90;
%  fi1 = 270;
%  fi2 = 273;%对结果没影响
%  fi3 = 270;
%  fi4 = 340;
%  fi5 = 330;
%  fi6 = 360;
%  k1 = 0.829;  %第一圆弧段造斜率
%  k2 = 2.007;  %第二圆弧段造斜率
%  k3 = 3.595;  %第三圆弧段造斜率
%  Dkop = 1000;  %第一造斜点垂深
%  Dd = 7000;  %第二造斜点垂深
%  Db = 10199.36;  %第三造斜点垂深
 
r1 = acos(cosd(a1) * cosd(a0) + sind(a1) * sind(a0) * cosd(fi2 - fi1));%第一造斜段（圆弧段）狗腿角
%r1 = 2 * asin(sqrt(sind( (a1- a0)/2 )^2 + sind( (fi2- fi1)/2 )^2 * sind(a0) * sind(a1) )); 

D1 = 100/(k1 * pi/180) * r1;%第一圆弧段长度

r2 = acos(cosd(a2) * cosd(a1) + sind(a2) * sind(a1) * cosd(fi4 - fi3));%第二造斜段（圆弧段）狗腿角
%r2 = 2 * asin(sqrt(sind( (a2- a1)/2 )^2 + sind( (fi4- fi3)/2 )^2 * sind(a1) * sind(a2) )); 

D3 = 100/(k2 * pi/180) * r2;%第二圆弧段长度

r3 = acos(cosd(a2) * cosd(a3) + sind(a2) * sind(a3) * cosd(fi6 - fi5));%第二造斜段（圆弧段）狗腿角
%r3 = 2 * asin(sqrt(sind( (a3- a2)/2 )^2 + sind( (fi6- fi5)/2 )^2 * sind(a2) * sind(a3) )); 

D5 = 100/(k3 * pi/180) * r3;%第三圆弧段长度
D2 = (Dd - Dkop - 100/(k1 * pi/180) * tan(r1/2) * (cosd(a1) + cosd(a0)))/cosd(a1);%第一稳斜段长度
D4 = (Db - Dd - 100/(k2 * pi/180) * tan(r2/2) * (cosd(a2) + cosd(a1)))/cosd(a2);%第一稳斜段长度
y=Dkop /cosd(a0) + D1 + D2 + D3 + D4 + D5 + 2500;%轨道总长度
end