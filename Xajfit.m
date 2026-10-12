function DC=Xajfit(x,allData,dataa)

date1 = datetime(1980, 1, 1, 8, 0, 0) ;
numSheets = length(allData);
DCZ=zeros(numSheets,1);

for i=1:numSheets

    dataZ=allData(i);
    dataZ=dataZ{1};
    data=dataZ(:,2:end);
    data=table2array(data);
    time=dataZ(:,1);
    time= table2array(time);
    time = datetime(time);
    %DCZ(j)=XXajfit(x,data);
    date2=time(1);
    time_diff = days(date2 - date1);
    days_ceil = floor(time_diff);
    index= days_ceil-90:days_ceil;
    datax=dataa(index,:);
    [Q,QZ]=XajfitQ(x,data,datax);
    NSE=jsdc(Q,QZ);
    SHICEHongfeng=max(QZ);
    MONIHongfeng=max(Q);
    SHICEFENGSHI=find(QZ==max(QZ));
    MONIFENGSHI=find(Q==max(Q));
    ERROROFFENGSHI=abs(MONIFENGSHI(1)-SHICEFENGSHI(1));  
    if ERROROFFENGSHI<=2
        wucha=0;
    elseif ERROROFFENGSHI<20
        wucha=(ERROROFFENGSHI-2)*0.1;
    else
        wucha=0;
    end
    DCZ(i)=1-NSE+abs((MONIHongfeng-SHICEHongfeng)/SHICEHongfeng)+wucha;
end

DC=sum(DCZ);
end