%Final 215 GPS Project
clc, clear
%Define GPS values
x1 = 16414028.668;
y1 = 660383.618; 
z1 = 20932036.907;
pr1 = 24658975.31743;

x2 = 16896800.648;
y2 = -18784061.365;
z2 = -7418318.856;
pr2 = 22964286.41228;

x3 = 9339639.616;
y3 = -14514964.658;
z3 = 20305107.161;
pr3 = 21338550.64536;

x4 = -18335582.591;
y4 = -11640868.305;
z4 = 15028599.071;
pr4 = 23606547.29359; 

%Compressed
a1 = 2*(x2 - x1); 
a2 = 2*(x3 - x2);
a3 = 2*(x4 - x3);

b1 = 2*(y2 - y1);
b2 = 2*(y3 - y2);
b3 = 2*(y4 - y3);

c1 = 2*(z2 - z1);
c2 = 2*(z3 - z2);
c3 = 2*(z4 - z3);

d1 = 2*(pr1 - pr2);
d2 = 2*(pr2 - pr3);
d3 = 2*(pr3 - pr4);

e1 = (x1^2 + y1^2 + z1^2) - (x2^2 + y2^2 + z2^2)- (pr1^2 - pr2^2); 
e2 = (x2^2 + y2^2 + z2^2) - (x3^2 + y3^2 + z3^2)- (pr2^2 - pr3^2);
e3 = (x3^2 + y3^2 + z3^2) - (x4^2 + y4^2 + z4^2)- (pr3^2 - pr4^2);


%Define Matrixes
M1 = [a1 b1 c1; a2 b2 c2; a3 b3 c3];
Md = [d1;d2;d3]; 
Me = [e1;e2;e3];

%Matrix inverse for division
M1inv = inv(M1);

MI = M1inv * -Md;  
MG = M1inv * -Me; 

%Final Compress

alpha2 = (MI(1,1)^2) + (MI(2,1)^2) + (MI(3,1)^2) - 1;
alpha1 = 2*(pr4 + MI(1,1)*(MG(1,1) - x4) + MI(2,1)*(MG(2,1) -y4) + MI(3,1)*(MG(3,1) -z4));
alpha0 = (MG(1,1) -x4)^2 + (MG(2,1) -y4)^2 + (MG(3,1) -z4)^2 - pr4^2;

temp = (alpha1^2 - 4*alpha2*alpha0);
Rc = (-alpha1 + sqrt( temp))./(2*alpha2);
%Rc2 = (-alpha1 - sqrt( temp))./(2*alpha2)
%x,y,z 
MI*Rc + MG
 




