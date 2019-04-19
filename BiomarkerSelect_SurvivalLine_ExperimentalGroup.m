%% signature for experiment group
load('~/TCGA_OV_DNAm.mat','DNAm_patient_name','ma','DNAm_pro');
load('~/TCGA_OV_clinic.mat','dmfs_day','dmfs_e','death_day','death_e','clinic_patient_name');
load('~/hub.mat');


[PatientName_intersect,ia,ib]=intersect(DNAm_patient_name,clinic_patient_name);
death_day_intersect=death_day(ib);
death_e_intersect=death_e(ib);
ma_intersect=ma(:,ia);
clear ia ib;


[C,ia,ib]=intersect(prebiomarker(:,1),DNAm_pro(:,1));
ma_prebiomarker=ma_intersect(ib,:);

clear ia C ma_intersect ib;

deathE_1=find(death_e_intersect==1);
deathE_1=deathE_1(randperm(length(deathE_1)));
deathE_1_HalfCount=length(deathE_1)/2;
train_group_deathE_1=deathE_1(1:deathE_1_HalfCount);
test_group_deathE_1=deathE_1(deathE_1_HalfCount+1:length(deathE_1));



deathE_0=find(death_e_intersect==0);
deathE_0=deathE_0(randperm(length(deathE_0)));
deathE_0_HalfCount=length(deathE_0)/2;
train_group_deathE_0=deathE_0(1:deathE_0_HalfCount);
test_group_deathE_0=deathE_0(deathE_0_HalfCount+1:length(deathE_0));



train_group=[train_group_deathE_1;train_group_deathE_0];
test_group=[test_group_deathE_1;test_group_deathE_0];

TrainGroup_PatientName=PatientName_intersect(train_group);
TestGroup_PatientName=PatientName_intersect(test_group);

save('C:\Users\���\Desktop\Test&Train_PatientName.mat','TrainGroup_PatientName','TestGroup_PatientName');

%clear dmfsE_1 dmfsE_1_HalfCount train_group_dmfsE_1 test_group_dmfsE_1 dmfsE_0 dmfsE_0_HalfCount train_group_dmfsE_0 test_group_dmfsE_0


%%  Cox

death_day_intersect_train=death_day_intersect(train_group);
death_e_intersect_train=death_e_intersect(train_group);
ma_train=ma_prebiomarker(:,train_group);
n=length(prebiomarker(:,1));
count=1;
biomarker_pos=[];
for i=1:n  
    data=(ma_train(i,:))';
    [~,~,~,stats_train]=coxphfit(data,death_day_intersect_train,'censoring',~death_e_intersect_train);
    %train
    if stats_train.p<0.05
    biomarker_pos=[biomarker_pos;i];
    result_train{count,1}=prebiomarker{i,1}; 
    result_train{count,2}=stats_train.beta; 
    result_train{count,3}=stats_train.p; %p
    count=count+1;
    end;
end;
count=count-1;



ma_SurviveLine_train=ma_train(biomarker_pos,:);
for j=1:count
    if result_train{j,2}<0
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

[con,ia,ib]=intersect(PatientName_intersect,expect_high_name);
one=ones(131);
one=one(1,:);
expect_high_final=[one;(death_day_intersect(ia))';(~death_e_intersect(ia))']';  

[con2,ia2,ib2]=intersect(PatientName_intersect,expect_low_name);
zero=zeros(132);
zero=zero(1,:);
expect_low_final=[zero;death_day_intersect(ia2)';(~death_e_intersect(ia2))']';

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

logrank_train(expect_final_train(expect_final_train(:,1)==0,2:3),expect_final_train(expect_final_train(:,1)==1,2:3)) 

clear hazard left right  stats


%%  selection
death_day_intersect_test=death_day_intersect(test_group);
death_e_intersect_test=death_e_intersect(test_group);
ma_SurviveLine_test=ma_prebiomarker(biomarker_pos,test_group);
for j=1:count
    if result_train{j,2}<0
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

[con,ia,ib]=intersect(PatientName_intersect,expect_high_name);
one=ones(131);
one=one(1,:);
expect_high_final=[one;(death_day_intersect(ia))';(~death_e_intersect(ia))']'; 

[con2,ia2,ib2]=intersect(PatientName_intersect,expect_low_name);
zero=zeros(132);
zero=zero(1,:);
expect_low_final=[zero;death_day_intersect(ia2)';(~death_e_intersect(ia2))']';

expect_final_test=[expect_low_final;expect_high_final];

clear colum_sum  patient_sum  expect_high_name  expect_low_name  one zero expect_high_final  expect_low_final con con2 ia ia2 ib ib2;


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

save('C:\Users\���\Desktop\hazard_result_train&test.mat','hazard_result_test','hazard_result_train');
logrank_test(expect_final_test(expect_final_test(:,1)==0,2:3),expect_final_test(expect_final_test(:,1)==1,2:3))  

clear hazard left right  stats



result_train=sortrows(result_train,3); 
xlswrite('C:\Users\���\Desktop\train_biomarker_final.xlsx',result_train);




















