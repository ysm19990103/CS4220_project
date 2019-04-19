%%  construction of dna methylation interaction network

load('~/TCGA_OV_DNAm.mat','DNAm_pro','ma');%%load data
ma_r=ma';
N=length(ma_r(1,:));
num=0;
for i=1:N  
    if sum(ma_r(:,i))~=0
        num=num+1;
    end;
end;


loc=(1:N)'; 
res=zeros(num*10,4);  

cou=0;
for i=1:N
    if sum(ma_r(:,i)~=0)  
        cou=cou+1;
        [co,p]=corr(ma_r,ma_r(:,i));
        da=[loc co p];
        da=da(~isnan(da(:,2)),:); 
        da=sortrows(da,-2); 
        res(((cou-1)*10+1):cou*10,1)=i;
        res(((cou-1)*10+1):cou*10,2:4)=da(2:11,1:3); 
    end;
end;

cou=1;
cou_max=length(res(:,1));
while cou<=cou_max
    result{cou,1}=DNAm_pro{res(cou,1),1};
    result{cou,2}=DNAm_pro{res(cou,2),1};
    result{cou,3}=res(cou,3);  
    result{cou,4}=res(cou,4);  
    cou=cou+1;
end

save('TCGA_DNAm_net_GBM.mat','result');






