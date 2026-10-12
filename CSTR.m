function CSZT=CSTR(x,E0,PALL)
 K=x(1);      
 WUM=x(2);     
 WLM=x(3);     
 WDM=x(4)-WLM-WUM;    
 C=x(5);     
 B=x(6);     
 IM=x(7);
 WM=WUM+WLM+WDM;
 EP=K*E0; 
 WMM=(1+B)*WM;
SZ=length(E0); 
FR=zeros(SZ,1);                        
EU=zeros(SZ,1);                        
EL=zeros(SZ,1);                       
ED=zeros(SZ,1);                       
E=zeros(SZ,1);                         
WU=zeros(SZ,1);                         
WL=zeros(SZ,1);                         
WD=zeros(SZ,1);                        
W=zeros(SZ,1);                        
PE=zeros(SZ,1); 
R=zeros(SZ,1);                  
WU(1)=WUM;                    
WL(1)=WLM;WD(1)=WDM;         
W(1)=WU(1)+WL(1)+WD(1);                
P=(1-IM)*PALL;                       
for I=1:SZ

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
a=WMM*(1-(1-W(I)/WM)^(1/(B+1)));  
if PE(I)<=0
    R(I)=0;
else
    if a+PE(I)<=WMM
    R(I)=PE(I)+W(I)-WM+WM*(1-(PE(I)+a)/WMM)^(B+1);
    else
    R(I)=PE(I)+W(I)-WM;
    end
end
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
end
CSZT=[WU(end),WL(end),WD(end)];

