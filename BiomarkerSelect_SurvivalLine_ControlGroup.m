load('~/TCGA_OV_DNAm.mat','DNAm_patient_name','ma','DNAm_pro');
load('~/TCGA_OV_clinic.mat','death_day','death_e','clinic_patient_name');
load('~/Test&Train_PatientName.mat');

%% control signature construction



[C_train,ia_train,ib_train]=intersect(TrainGroup_PatientName,clinic_patient_name);
death_day_intersect_train=death_day(ib_train);
death_e_intersect_train=death_e(ib_train);

[C,ia,ib]=intersect(TrainGroup_PatientName,DNAm_patient_name);
ma_train=ma(:,ib);

clear ia ib;






[C_test,ia_test,ib_test]=intersect(TestGroup_PatientName,clinic_patient_name);
death_day_intersect_test=death_day(ib_test);
death_e_intersect_test=death_e(ib_test);

[C,ia,ib]=intersect(TestGroup_PatientName,DNAm_patient_name);
ma_test=ma(:,ib);

clear ia ib;





%Cox

n=length(DNAm_pro(:,1));
for i=1:n  
    data=(ma_train(i,:))';
    [~,~,~,stats_train]=coxphfit(data,death_day_intersect_train,'censoring',~death_e_intersect_train);%��i������ı��ֵ�����֮����cox�ع�
   
    
    result_train{i,1}=i; 
    result_train{i,2}=DNAm_pro{i,1}; 
    result_train{i,3}=stats_train.beta; 
    result_train{i,4}=stats_train.p; %p-value
 end;
 
 result_train=sortrows(result_train,4);
 result_train=result_train(1:76,:);
 biomarker_pos=result_train(1:76,1);
 
 
 
 



%% selection
biomarker_pos=cell2mat(biomarker_pos);
ma_SurviveLine_train=ma_train(biomarker_pos,:);
for j=1:length(biomarker_pos)
    if result_train{j,3}<0
        ma_SurviveLine_train(j,:)=(-1)*(ma_SurviveLine_train(j,:));
    end
   %ma_SurviveLine_train(j,:)=result_train{j,2}.*(ma_SurviveLine_train(j,:));
end

colum_sum=sum(ma_SurviveLine_train);
for i=1:length(colum_sum)
   patient_sum{i,1}=TrainGroup_PatientName{i,1};
   patient_sum{i,2}=colum_sum(1,i);
end
 result=sortrows(patient_sum,2);

expect_low_name=result(1:ceil(length(colum_sum)/2),1);
expect_high_name=result((ceil(length(colum_sum)/2)+1):length(colum_sum),1);

[con,ia,ib]=intersect(TrainGroup_PatientName,expect_high_name);
one=ones(131);
one=one(1,:);
expect_high_final=[one;(death_day_intersect_train(ia))';(~death_e_intersect_train(ia))']';  

[con2,ia2,ib2]=intersect(TrainGroup_PatientName,expect_low_name);
zero=zeros(132);
zero=zero(1,:);
expect_low_final=[zero;death_day_intersect_train(ia2)';(~death_e_intersect_train(ia2))']';

expect_final_train=[expect_low_final;expect_high_final];


clear colum_sum  patient_sum  expect_high_name  expect_low_name  one zero expect_high_final  expect_low_final con con2 ia ia2 ib ib2;


[~,~,~,stats]=coxphfit(expect_final_train(:,1),expect_final_train(:,2),'censoring',expect_final_train(:,3));
hazard=exp(stats.beta);              %Hazard
left=exp(stats.beta-1.96*stats.se);  %0.95CI left
right=exp(stats.beta+1.96*stats.se); %0.95 CI right

hazard_result_train{1,1}=hazard;  
hazard_result_train{1,2}=left;     
hazard_result_train{1,3}=right;
hazard_result_train{1,4}=logrank_p(expect_final_train(expect_final_train(:,1)==0,2:3),expect_final_train(expect_final_train(:,1)==1,2:3));  %p-value
hazard_result_train{1,5}=length(expect_final_train(:,1))-sum(expect_final_train(:,1));  
hazard_result_train{1,6}=sum(expect_final_train(:,1));                    

logrank_control_train(expect_final_train(expect_final_train(:,1)==0,2:3),expect_final_train(expect_final_train(:,1)==1,2:3))  

clear hazard left right  stats



ma_SurviveLine_test=ma_test(biomarker_pos,:);
for j=1:length(biomarker_pos)
    if result_train{j,3}<0
        ma_SurviveLine_test(j,:)=(-1)*(ma_SurviveLine_test(j,:));
    end
     %ma_SurviveLine_test(j,:)=result_train{j,2}.*(ma_SurviveLine_test(j,:));
end
colum_sum=sum(ma_SurviveLine_test);
for i=1:length(colum_sum)
   patient_sum{i,1}=TestGroup_PatientName{i,1};
   patient_sum{i,2}=colum_sum(1,i);
end
result=sortrows(patient_sum,2);

expect_low_name=result(1:ceil(length(colum_sum)/2),1);
expect_high_name=result((ceil(length(colum_sum)/2)+1):length(colum_sum),1);

[con,ia,ib]=intersect(TestGroup_PatientName,expect_high_name);
one=ones(131);
one=one(1,:);
expect_high_final=[one;(death_day_intersect_test(ia))';(~death_e_intersect_test(ia))']'; 

[con2,ia2,ib2]=intersect(TestGroup_PatientName,expect_low_name);
zero=zeros(132);
zero=zero(1,:);
expect_low_final=[zero;death_day_intersect_test(ia2)';(~death_e_intersect_test(ia2))']';

expect_final_test=[expect_low_final;expect_high_final];

clear colum_sum  patient_sum  expect_high_name  expect_low_name  one zero expect_high_final  expect_low_final con con2 ia ia2 ib ib2 i j ;



[~,~,~,stats]=coxphfit(expect_final_test(:,1),expect_final_test(:,2),'censoring',expect_final_test(:,3));
hazard=exp(stats.beta);              %Hazard
left=exp(stats.beta-1.96*stats.se);  %0.95CI left
right=exp(stats.beta+1.96*stats.se); %0.95 CI right

hazard_result_test{1,1}=hazard;  
hazard_result_test{1,2}=left;     
hazard_result_test{1,3}=right;
hazard_result_test{1,4}=logrank_p(expect_final_test(expect_final_test(:,1)==0,2:3),expect_final_test(expect_final_test(:,1)==1,2:3));  %p-value
hazard_result_test{1,5}=length(expect_final_test(:,1))-sum(expect_final_test(:,1)); 
hazard_result_test{1,6}=sum(expect_final_test(:,1));                   

save('~/hazard_result_train&test.mat','hazard_result_test','hazard_result_train');
logrank_control_test(expect_final_test(expect_final_test(:,1)==0,2:3),expect_final_test(expect_final_test(:,1)==1,2:3)) 


