function Q=MSJG(Qin,xx,Km,t)
SZ=length(Qin);
Q=zeros(SZ,1);
C0 = (0.5*t - Km*xx)/(0.5*t + Km - Km*xx);
C1 = (0.5*t + Km*xx)/(0.5*t + Km - Km*xx);
C2 = (-0.5*t + Km - Km*xx)/(0.5*t + Km - Km*xx);


Q(1)=Qin(1);
 for I=2:SZ
    Q(I)=C0*Qin(I)+C1*Qin(I-1)+C2*Q(I-1);
 end
end
