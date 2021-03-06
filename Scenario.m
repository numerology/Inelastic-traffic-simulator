function [bs] = Scenario(intersiteD)

x_hexagon=[-1 -0.5 0.5 1 0.5 -0.5 -1];
y_hexagon=[0 -sqrt(3)/2 -sqrt(3)/2 0 sqrt(3)/2 sqrt(3)/2 0];
multiple=intersiteD/(2*cosd(30));
s1=[1,0];
s2=[-0.5,sqrt(3)/2];
s3=[-0.5,-sqrt(3)/2];


%% bs mat
bs(1,:)=[(multiple*(+3*0)),multiple*(+sqrt(3)*2)];%1
bs(2,:)=[(multiple*(-1.5+3*0)),multiple*(+sqrt(3)/2+sqrt(3)*1)];%2
bs(3,:)=[(multiple*(+3*0)),multiple*(+sqrt(3)*1)];%3
bs(4,:)=[(multiple*(+1.5+3*0)),multiple*(+sqrt(3)/2+sqrt(3)*1)];%4
bs(5,:)=[(multiple*(+3*-1)),multiple*(+sqrt(3)*1)];%5
bs(6,:)=[(multiple*(-1.5+3*0)),multiple*(-sqrt(3)/2+sqrt(3)*1)];%6
bs(7,:)=[(multiple*(+3*0)),multiple*(+sqrt(3)*0)];%7
bs(8,:)=[(multiple*(+1.5+3*0)),multiple*(+sqrt(3)/2+sqrt(3)*0)];%8
bs(9,:)=[(multiple*(+3*1)),multiple*(+sqrt(3)*1)];%9
bs(10,:)=[(multiple*(+3*-1)),multiple*(+sqrt(3)*0)];%10
bs(11,:)=[(multiple*(-1.5+3*0)),multiple*(-sqrt(3)/2+sqrt(3)*0)];%11
bs(12,:)=[(multiple*(+3*0)),multiple*(+sqrt(3)*-1)];%12
bs(13,:)=[(multiple*(+1.5+3*0)),multiple*(+sqrt(3)/2+sqrt(3)*-1)];%13
bs(14,:)=[(multiple*(+3*1)),multiple*(+sqrt(3)*0)];%14
bs(15,:)=[(multiple*(+3*-1)),multiple*(+sqrt(3)*-1)];%15
bs(16,:)=[(multiple*(-1.5+3*0)),multiple*(-sqrt(3)/2+sqrt(3)*-1)];%16
bs(17,:)=[(multiple*(+3*0)),multiple*(+sqrt(3)*-2)];%17
bs(18,:)=[(multiple*(+1.5+3*0)),multiple*(+sqrt(3)/2+sqrt(3)*-2)];%18
bs(19,:)=[(multiple*(+3*1)),multiple*(+sqrt(3)*-1)];%19

bs(20,:)=[(multiple*(+3*0)),multiple*(+sqrt(3)*3)];%20
bs(21,:)=[(multiple*(+1.5+3*0)),multiple*(+sqrt(3)/2+sqrt(3)*2)];%21
bs(22,:)=[(multiple*(+3*1)),multiple*(+sqrt(3)*2)];%22
bs(23,:)=[(multiple*(+1.5+3*1)),multiple*(+sqrt(3)/2+sqrt(3)*1)];%23
bs(24,:)=[(multiple*(-1.5+3*2)),multiple*(-sqrt(3)/2+sqrt(3)*1)];%24
bs(25,:)=[(multiple*(+1.5+3*1)),multiple*(+sqrt(3)/2+sqrt(3)*-1)];%25
bs(26,:)=[(multiple*(+1.5+3*1)),multiple*(+sqrt(3)/2+sqrt(3)*-2)];%26
bs(27,:)=[(multiple*(+3*1)),multiple*(+sqrt(3)*-2)];%27
bs(28,:)=[(multiple*(+1.5+3*0)),multiple*(+sqrt(3)/2+sqrt(3)*-3)];%28
bs(29,:)=[(multiple*(+3*0)),multiple*(+sqrt(3)*-3)];%29
bs(30,:)=[(multiple*(-1.5+3*0)),multiple*(-sqrt(3)/2+sqrt(3)*-2)];%30
bs(31,:)=[(multiple*(+3*-1)),multiple*(+sqrt(3)*-2)];%31
bs(32,:)=[(multiple*(+1.5+3*-2)),multiple*(+sqrt(3)/2+sqrt(3)*-2)];%32
bs(33,:)=[(multiple*(+1.5+3*-2)),multiple*(+sqrt(3)/2+sqrt(3)*-1)];%33
bs(34,:)=[(multiple*(+1.5+3*-2)),multiple*(+sqrt(3)/2+sqrt(3)*0)];%34
bs(35,:)=[(multiple*(+1.5+3*-2)),multiple*(+sqrt(3)/2+sqrt(3)*1)];%35
bs(36,:)=[(multiple*(+3*-1)),multiple*(+sqrt(3)*2)];%36
bs(37,:)=[(multiple*(-1.5+3*0)),multiple*(+sqrt(3)/2+sqrt(3)*2)];%38

bs(38,:)=[multiple*(+3*0),multiple*(+sqrt(3)*4)];%38
bs(39,:)=[multiple*(+1.5+3*0),multiple*(+sqrt(3)/2+sqrt(3)*3)];%39
bs(40,:)=[multiple*(+3*1),multiple*(+sqrt(3)*3)];%40
bs(41,:)=[multiple*(+1.5+3*1),multiple*(+sqrt(3)/2+sqrt(3)*2)];%41
bs(42,:)=[multiple*(+3*2),multiple*(+sqrt(3)*2)];%42
bs(43,:)=[multiple*(+3*2),multiple*(+sqrt(3)*1)];%43
bs(44,:)=[multiple*(+3*2),multiple*(+sqrt(3)*0)];%44
bs(45,:)=[multiple*(+3*2),multiple*(+sqrt(3)*-1)];%45
bs(46,:)=[multiple*(+3*2),multiple*(+sqrt(3)*-2)];%46
bs(47,:)=[multiple*(+1.5+3*1),multiple*(+sqrt(3)/2+sqrt(3)*-3)];%47
bs(48,:)=[multiple*(+3*1),multiple*(+sqrt(3)*-3)];%48
bs(49,:)=[multiple*(+1.5+3*0),multiple*(+sqrt(3)/2+sqrt(3)*-4)];%49
bs(50,:)=[multiple*(+3*0),multiple*(+sqrt(3)*-4)];%50
bs(51,:)=[multiple*(-1.5+3*0),multiple*(+sqrt(3)/2+sqrt(3)*-4)];%51
bs(52,:)=[multiple*(+3*-1),multiple*(+sqrt(3)*-3)];%52
bs(53,:)=[multiple*(+1.5+3*-2),multiple*(+sqrt(3)/2+sqrt(3)*-3)];%53
bs(54,:)=[multiple*(+3*-2),multiple*(+sqrt(3)*-2)];%54
bs(55,:)=[multiple*(+3*-2),multiple*(+sqrt(3)*-1)];%55
bs(56,:)=[multiple*(+3*-2),multiple*(+sqrt(3)*0)];%56
bs(57,:)=[multiple*(+3*-2),multiple*(+sqrt(3)*1)];%57
bs(58,:)=[multiple*(+3*-2),multiple*(+sqrt(3)*2)];%58
bs(59,:)=[multiple*(+1.5+3*-2),multiple*(+sqrt(3)/2+sqrt(3)*2)];%59
bs(60,:)=[multiple*(+3*-1),multiple*(+sqrt(3)*3)];%60
bs(61,:)=[multiple*(-1.5+3*0),multiple*(+sqrt(3)/2+sqrt(3)*3)];%61

bsextra=[1:19,15,16,17,5,10,15,1,2,5,9,4,1,19,14,9,17,18,19,10,11,12,18,2,6,11,16,4,3,6,10,14,8,3,2,18,13,8,4,16,12,13,14];
bs=[bs,bsextra'];

