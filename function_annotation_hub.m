   load('~/TCGA_OV_DNAm.mat','probe');
   result={};
   select_feature_A=hub_probe;
   for j=1:length(KEGG(:,1))
        count=0;
        for k=2:length(KEGG(1,:))
            if(~isequal(KEGG{j,k},''))
                count=count+1;
            end
        end
        select_feature_B=KEGG(j,2:count);
        select_feature_C=intersect(select_feature_A,select_feature_B);
        p_value= 1-hygecdf(length(select_feature_C)-1,length(probe),length(select_feature_A),length(select_feature_B));
        
        
        result{j,1}=KEGG{j,1};
        result{j,2}=length(select_feature_C);
        result{j,3}=p_value;
        
        clear p_value;
        clear select_feature_B;
        clear select_feature_C;
   end
   result=sortrows(result,3);