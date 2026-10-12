function Q=Ziliuyu(x,E0,F,PALL,Q0,cszt)
 K=x(1);      
 WUM=x(2);            
 WLM=x(3);     
 WDM=x(4)-WUM-WLM;       
 C=x(5);       
 B=x(6);     
 IM=x(7);  
 SM=x(8);
 EX = x(9); 
 KG=x(10); 
 KI = x(11);
 axx=x(12);
 bxx=x(13);
c_min=min(x(14),x(15));
c_max=max(x(14),x(15));

aIxx=x(16);
bIxx=x(17);
cI_min=min(x(18),x(19));
cI_max=max(x(18),x(19));

aGxx=x(20);
bGxx=x(21);
cG_min=min(x(22),x(23));
cG_max=max(x(22),x(23));

      
  

 t = 1;                  


 U=F/(3.6*t);
 WM=WUM+WLM+WDM;
 EP=K*E0; 
 WMM=(1+B)*WM;
SMM=(1+EX)*SM;
SZ=length(E0); 
S1=zeros(SZ,1);S1(1)=0;               %为S1预分配内存；S1表示本时段初产流面积上的平均自由水深
FR=zeros(SZ,1);                        %为FR预分配内存；FR表示本时段产流面积比例
EU=zeros(SZ,1);                        %为EU预分配内存
EL=zeros(SZ,1);                        %为EL预分配内存
ED=zeros(SZ,1);                        %为ED预分配内存
E=zeros(SZ,1);                         %为E预分配内存
WU=zeros(SZ,1);                         %为WU预分配内存
WL=zeros(SZ,1);                         %为WL预分配内存
WD=zeros(SZ,1);                         %为WD预分配内存
W=zeros(SZ,1);                         %为W预分配内存
R=zeros(SZ,1);                         %为R预分配内存
RS3=zeros(SZ,1);                         %为RS3预分配内存；RS3表示三水源地面径流
RI3=zeros(SZ,1);                         %为RI3预分配内存；RI3表示三水源壤中水径流
RG3=zeros(SZ,1);                         %为RG3预分配内存；RG3表示三水源地下径流
PE=zeros(SZ,1);                         %为PE预分配内存；PE表示净雨量
QS3=zeros(SZ,1);  
QI3=zeros(SZ,1);                       %为QI3预分配内存；QI3表示三分水源壤中流汇流
QG3=zeros(SZ,1);                        %为QG3预分配内存；QG3表示三分水源地下径流汇流

  
                    %为Qz预分配内存；Qz表示各个子流域出口断面流量过程
WU(1)= cszt(1);                     %流域上层初始土壤含水量   范围：0～20
WL(1)= cszt(2);
WD(1)= cszt(3);         %WL、WD都表示本时段初相应土层土壤含水量
W(1)=WU(1)+WL(1)+WD(1);                %W表示本时段初土壤含水量
QG3(1)=Q0;
CSA=WMM*(1-(1-W(1)/WM)^(1/(B+1))); 
FR0=1-(1-CSA/WMM)^B;

P=(1-IM)*PALL;   %P表示某子流域降雨过程，data数组中第4列到第10列分别为
                        % 子流域1-7的降雨资料
RB=IM*PALL;                %RB表示降落在不透水面上的降雨量
for I=1:SZ
%%%%%%%%%%%%%三层蒸散发模型%%%%%%%%%%%%%
if WU(I)+P(I)>=EP(I)
    EU(I)=EP(I);EL(I)=0;ED(I)=0;
else
    if WL(I)>=C*WLM
        EU(I)=WU(I)+P(I);EL(I)=(EP(I)-EU(I))*WL(I)/WLM;ED(I)=0;
    else
        if WL(I)>=C*(EP(I)-(WU(I)+P(I)))
            EU(I)=WU(I)+P(I);EL(I)=C*(EP(I)-EU(I));ED(I)=0;
        else
            EU(I)=WU(I)+P(I);EL(I)=WL(I);ED(I)=C*(EP(I)-EU(I))-EL(I);
        end
    end
end
E(I)=EU(I)+EL(I)+ED(I);PE(I)=P(I)-E(I);
%%%%%%%%%%%%%三层蒸散发模型%%%%%%%%%%%%%
%%%%%%%%%%%%%计算总径流深%%%%%%%%%%%%%
a=WMM*(1-(1-W(I)/WM)^(1/(B+1)));  %a表示流域初始土壤含水量最大值
if PE(I)<=0
    R(I)=0;
else
    if a+PE(I)<=WMM
    R(I)=PE(I)+W(I)-WM+WM*(1-(PE(I)+a)/WMM)^(B+1);
    else
    R(I)=PE(I)+W(I)-WM;
    end
end
%%%%%%%%%%%%%计算总径流深%%%%%%%%%%%%%
%%%%%%%%%%%%%计算下一时段土壤含水量%%%%%%%%%%%%%
WU(I+1)=WU(I)+P(I)-EU(I)-R(I);
WL(I+1)=WL(I)-EL(I);
WD(I+1)=WD(I)-ED(I);
if WU(I+1)>WUM
    WL(I+1)=WL(I)-EL(I)+WU(I+1)-WUM;
    WU(I+1)=WUM;
end
if WL(I+1)>WLM
    WD(I+1)=WD(I)-ED(I)+WL(I+1)-WLM;
    WL(I+1)=WLM;
end
if WD(I+1)>WDM
    WD(I+1)=WDM;
end
W(I+1)=WU(I+1)+WL(I+1)+WD(I+1);
LH=WMM*(1-(1-W(I+1)/WM)^(1/(B+1)));
FR(I)=1-(1-LH/WMM)^B;
%%%%%%%%%%%%%计算下一时段土壤含水量%%%%%%%%%%%%%
%%%%%%%%%%%%%计算产流面积FR=R/PE%%%%%%%%%%%%%

%%%%%%%%%%%%%计算产流面积FR=R/PE%%%%%%%%%%%%%
%%%%%%%%%%%%%三分水源计算%%%%%%%%%%%%%
 %%%%%%%%%%%%%第一时段数据计算%%%%%%%%%%%%%
 if I==1
      if R(I)>0
          if (S1(I)*FR0/FR(I))/SM>1
              AU = SMM;
          else
          AU=SMM*(1-(1-(S1(I)*FR0/FR(I))/SM)^(1/(1+EX)));
          end
          if R(I)/FR(I)+AU<SMM
              RS3(I)=FR(I)*(R(I)/FR(I)+S1(I)*FR0/FR(I)-SM+SM*(1-(R(I)/FR(I)+AU)/SMM)^(EX+1));
              S=S1(I)*FR0/FR(I)+(R(I)-RS3(I))/FR(I);
              RI3(I)=KI*S*FR(I);
              RG3(I)=KG*S*FR(I);
              S1(I+1)=S*(1-KI-KG);
          else
               RS3(I)=FR(I)*(R(I)/FR(I)+S1(I)*FR0/FR(I)-SM );
               S=S1(I)*FR0/FR(I)+(R(I)-RS3(I))/FR(I);
               RI3(I)=KI*S*FR(I);
               RG3(I)=KG*S*FR(I);
               S1(I+1)=S*(1-KI-KG);
          end
      else
          S=S1(I)*FR0/FR(I);
          RS3(I)=0;RI3(I)=KI*S*FR(I);RG3(I)=KG*S*FR(I);S1(I+1)=S*(1-KI-KG);  
      end
%%%%%%%%%%%%%第一时段数据计算%%%%%%%%%%%%%
%%%%%%%%%%%%%其余时段数据计算%%%%%%%%%%%%%
 else
     if R(I)>0
          if (S1(I)*FR(I-1)/FR(I))/SM>1
              AU = SMM;
          else
              AU=SMM*(1-(1-(S1(I)*FR(I-1)/FR(I))/SM)^(1/(1+EX)));
          end
          if R(I)/FR(I)+AU<SMM
              RS3(I)=FR(I)*(R(I)/FR(I)+S1(I)*FR(I-1)/FR(I)-SM+SM*(1-(R(I)/FR(I)+AU)/SMM)^(EX+1));
              S=S1(I)*FR(I-1)/FR(I)+(R(I)-RS3(I))/FR(I);
              RI3(I)=KI*S*FR(I);
              RG3(I)=KG*S*FR(I);
              S1(I+1)=S*(1-KI-KG);
          else
              RS3(I)=FR(I)*(R(I)/FR(I)+S1(I)*FR(I-1)/FR(I)-SM);
              S=S1(I)*FR(I-1)/FR(I)+(R(I)-RS3(I))/FR(I);
              RI3(I)=KI*S*FR(I);
              RG3(I)=KG*S*FR(I);
              S1(I+1)=S*(1-KI-KG);
          end
     else
            S=S1(I)*FR(I-1)/FR(I);
            RS3(I)=0;RI3(I)=KI*S*FR(I);RG3(I)=KG*S*FR(I);S1(I+1)=S*(1-KI-KG);  
     end
 end
              
        
     
        
end

%%%%%%%%%%%%%其余时段数据计算%%%%%%%%%%%%%
%%%%%%%%%%%%%三分水源计算%%%%%%%%%%%%%

%%%%%%%%%%%%%三分水源汇流计算%%%%%%%%%%%%%
QSS3=(RS3+RB);                               %QS3为地面径流的坡地汇流

QS3= dynamic_reservoir_currentS_Q0isQS1(QSS3,U,t,QS3(1),c_min,c_max,axx,bxx);
QI3= dynamic_reservoir_currentS_Q0isQS1(RI3,U,t,QI3(1),cI_min,cI_max,aIxx,bIxx);
QG3= dynamic_reservoir_currentS_Q0isQS1(RG3,U,t,QG3(1),cG_min,cG_max,aGxx,bGxx);




Q=QS3+QI3+QG3;


