%load('~/BP.mat');
%load('~/DNAm_pro_selected_finally.mat');
%load('~/row_num_later.mat');
%load('~/TCGA_OV_DNAm.mat','probe');

result_final={};
for i=1:length(DNAm_pro_selected(:,1))
    
    select_feature_A=DNAm_pro_selected(i,1:row_num(i,1));
    
    
    for j=1:length(BP(:,1))
        count=0;
        for k=2:length(BP(1,:))
            if(~isequal(BP{j,k},''))
                count=count+1;
            end
        end
        select_feature_B=BP(j,2:count);
        select_feature_C=intersect(select_feature_A,select_feature_B);
        p_value= 1-hygecdf(length(select_feature_C)-1,length(probe),length(select_feature_A),length(select_feature_B));
        
        result{j,1}=i;
        result{j,2}=BP{j,1};
        result{j,3}=length(select_feature_C);
        result{j,4}=p_value;
        
        clear p_value;
        clear select_feature_B;
        clear select_feature_C;
    end
        
    result=sortrows(result,4);
    result=result(1:10,1:4);
    
    begin=(i-1)*10+1;
    ending=(i-1)*10+10;
    
    result_final(begin:ending,1:4)=result;
   
       
    clear select_feature_A;
    clear result;
    
end
    
xlswrite('~/final_result.xlsx',result_final);