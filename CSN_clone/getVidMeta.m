function [HS_Data] = getVidMeta(folder)


fn = getfn(folder,'config.xsv$');

for i = 1:length(fn)
    filename=fn(1,i);
    fileID=fopen(filename{1,1});
    data = textscan(fileID,'%s',40,'Delimiter','\n');
    data=data{1,1};
    data=data(2:end);
    for j = 1:length(data)
        tosplit=data(j);
        C = strsplit(tosplit{1,1},'=');
        FieldName = strrep(cell2mat(C(1)),' ','_');
        FieldName = strrep(FieldName,'-','_');
        Metadata.(FieldName)=C(2);
    end
    HS_Metadata{1,i}=Metadata;
    fclose(fileID);
end


%%

VidNames=[];
Hour=[];
Minute=[];
Second=[];
Milliseconds=[];
Time=[];
NFrames=[];
MinIdx=[];
StartIdx=[];
StopIdx=[];
Rate=[];
Date=[];
for vid = 1:length(HS_Metadata)
    VidNames= [VidNames; HS_Metadata{1,vid}.Name{1,1}];
    NFrames= [NFrames; str2num(HS_Metadata{1,vid}.Frames{1,1})];
    MinIdx=[MinIdx;str2num(HS_Metadata{1,vid}.MinIdx{1,1})];
    StartIdx=[StartIdx;str2num(HS_Metadata{1,vid}.StartIdx{1,1})];
    StopIdx=[StopIdx;str2num(HS_Metadata{1,vid}.StopIdx{1,1})];
    Rate=[Rate;str2num(HS_Metadata{1,vid}.Rate{1,1})];
    d=HS_Metadata{1,vid}.Day{1,1};
    if length(d==1)
        d=strcat('0',d);
    end
    mo=HS_Metadata{1,vid}.Month{1,1};
    if length(mo)==1
        mo=strcat('0',mo);
    end
    y=HS_Metadata{1,vid}.Year{1,1};
    Date=[Date; strcat(d,'-',mo,'-',y)];
    Hour=[Hour;str2num(HS_Metadata{1,vid}.Hour{1,1})];
    Minute=[Minute;str2num(HS_Metadata{1,vid}.Minute{1,1}) ];
    Second=[Second;str2num(HS_Metadata{1,vid}.Second{1,1}) ];
    Milliseconds=[Milliseconds;str2num(HS_Metadata{1,vid}.Milliseconds{1,1}) ];
    h=HS_Metadata{1,vid}.Hour{1,1};
    if length(h)==1
        h=strcat('0',h);
    end
    m=HS_Metadata{1,vid}.Minute{1,1};
    if length(m)==1
        m=strcat('0',m);
    end
    s=HS_Metadata{1,vid}.Second{1,1};
    if length(s)==1
        s=strcat('0',s);
    end
    ms=HS_Metadata{1,vid}.Milliseconds{1,1};
    if length(ms)==1
        ms=strcat('0',ms);
    end
    if length(ms)==2
        ms=strcat('0',ms);
    end
    t=strcat(h,':',m,':',s,'.',ms);
    Time=[Time;t];

end

HS_Data.VidNames=VidNames;
HS_Data.Time=Time;
HS_Data.NFrames=NFrames;
HS_Data.MinIdx=MinIdx;
HS_Data.StartIdx=StartIdx;
HS_Data.StopIdx=StopIdx;
HS_Data.Rate=Rate;
HS_Data.Date=Date;

