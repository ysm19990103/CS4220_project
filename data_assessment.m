
%% data stability assessment

p_value=[p_value_DNAm;p_value_mRNA;p_value_miRNA];
for i=1:100
    group{i,1}='DNA methylation';
end;
for i=101:200
    group{i,1}='mRNA expression';
end;
for i=201:300
    group{i,1}='miRNA expression'; 
end;
h=boxplot(p_value,group);
p1=ranksum(p_value_DNAm,p_value_mRNA);  %
p2=ranksum(p_value_DNAm,p_value_miRNA); %
save('~/p_final(2).mat','p1','p2','p_value_DNAm','p_value_mRNA','p_value_miRNA')
saveas(h,'~/finalresult(2).jpg') 














