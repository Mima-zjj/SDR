function DC=jsdc(Q,QZ)
d=0;dz=0;
SZ=length(Q); 
for I=1:SZ
    d=d+Q(I);
    dz=dz+QZ(I);
end
D=100*(d-dz)/dz;            %计算年径流深相对误差
RZ=dz*24*3600/(290.1*1000); %RZ表示实测年径流深
Rd=d*24*3600/(290.1*1000);  %Rd表示计算年径流深
c=0;dc=0;
for I=1:SZ
    c=c+(Q(I)-QZ(I))^2;
    dc=dc+(QZ(I)-dz/SZ)^2;
end
DC=1-c/dc;                 %计算确定性系数