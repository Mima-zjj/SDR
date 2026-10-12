function [Q,QZ]=XajfitQ(x,data,data1)
QZ=data(:,1);                        %QZ表示实测流量
E0=data(:,2);  
F=[70 74.4 60 101 189 248 114 36 122.8 268 169 159 176.5 26.3];
%PALL=data(:,3);
PALL1=data(:,4);
PALL2=data(:,5);
PALL3=data(:,6);
PALL4=data(:,7);
PALL5=data(:,8);
PALL6=data(:,9);
PALL7=data(:,10);
PALL8=data(:,11);
PALL9=data(:,12);
PALL10=data(:,13);
PALL11=data(:,14);
PALL12=data(:,15);
PALL13=data(:,16);
PALL14=data(:,17);

E0r=data1(:,2);
PALLr1=data1(:,4);
PALLr2=data1(:,5);
PALLr3=data1(:,6);
PALLr4=data1(:,7);
PALLr5=data1(:,8);
PALLr6=data1(:,9);
PALLr7=data1(:,10);
PALLr8=data1(:,11);
PALLr9=data1(:,12);
PALLr10=data1(:,13);
PALLr11=data1(:,14);
PALLr12=data1(:,15);
PALLr13=data1(:,16);

t=1;
%SZ=length(QZ);
%Q=zeros(SZ,1);
Q_chushi=QZ(1);
V=x(24);
xx1=x(25);
xx2=x(26);

 Q001=F(1)/sum(F)*Q_chushi;
 cstr1=CSTR(x,E0r,PALLr1);
 Q1=Ziliuyu(x,E0,F(1),PALL1,Q001,cstr1);
 
 Q002=F(2)/sum(F)*Q_chushi;
 cstr2=CSTR(x,E0r,PALLr2);
 Q2=Ziliuyu(x,E0,F(2),PALL2,Q002,cstr2); 
 
 Q003=F(3)/sum(F)*Q_chushi;
 cstr3=CSTR(x,E0r,PALLr3);
 Q3=Ziliuyu(x,E0,F(3),PALL3,Q003,cstr3);  
 
 Q004=F(4)/sum(F)*Q_chushi;
 cstr4=CSTR(x,E0r,PALLr4);
 Q4=Ziliuyu(x,E0,F(4),PALL4,Q004,cstr4);  
 
 Q005=F(5)/sum(F)*Q_chushi;
 cstr5=CSTR(x,E0r,PALLr5);
 Q5=Ziliuyu(x,E0,F(5),PALL5,Q005,cstr5);  
 
 Q006=F(6)/sum(F)*Q_chushi;
 cstr6=CSTR(x,E0r,PALLr6);
 Q6=Ziliuyu(x,E0,F(6),PALL6,Q006,cstr6);  
 
 Q007=F(7)/sum(F)*Q_chushi;
 cstr7=CSTR(x,E0r,PALLr7);
 Q7=Ziliuyu(x,E0,F(7),PALL7,Q007,cstr7);  
 
 Q008=F(8)/sum(F)*Q_chushi;
 cstr8=CSTR(x,E0r,PALLr8);
 Q8=Ziliuyu(x,E0,F(8),PALL8,Q008,cstr8); 

 Q009=F(9)/sum(F)*Q_chushi;
 cstr9=CSTR(x,E0r,PALLr9);
 Q9=Ziliuyu(x,E0,F(9),PALL9,Q009,cstr9); 

 Q010=F(10)/sum(F)*Q_chushi;
 cstr10=CSTR(x,E0r,PALLr10);
 Q10=Ziliuyu(x,E0,F(10),PALL10,Q010,cstr10); 

 Q011=F(11)/sum(F)*Q_chushi;
 cstr11=CSTR(x,E0r,PALLr11);
 Q11=Ziliuyu(x,E0,F(11),PALL11,Q011,cstr11);

 Q012=F(12)/sum(F)*Q_chushi;
 cstr12=CSTR(x,E0r,PALLr12);
 Q12=Ziliuyu(x,E0,F(12),PALL12,Q012,cstr12);

 Q013=F(13)/sum(F)*Q_chushi;
 cstr13=CSTR(x,E0r,PALLr13);
 Q13=Ziliuyu(x,E0,F(13),PALL13,Q013,cstr13);  
 
 Q14=PALL14*F(14)/(3.6*t);
Q14(1)=F(14)/sum(F)*Q_chushi;
KmZ36=41495/V/3600;
KmZ45=25029/V/3600;

 N_36=round(KmZ36);

 if N_36==0
     disp("子流域36汇流时段过短")
 else
     Km36=KmZ36/N_36;
     for i36=1:N_36
         Q3=MSJG(Q3,xx1,Km36,t);
         Q6=MSJG(Q6,xx1,Km36,t);
     end
 end

 N_45=round(KmZ45);

 if N_45==0
     disp("子流域45汇流时段过短")
 else
     Km45=KmZ45/N_45;
     for i45=1:N_45
         Q4=MSJG(Q4,xx2,Km45,t);
         Q5=MSJG(Q5,xx2,Km45,t);
     end
 end

Q=Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8+Q9+Q10+Q11+Q12+Q13+Q14;


end

   




