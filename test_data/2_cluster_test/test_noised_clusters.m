nsamples = 100;
step = 1;
noisepoints = 100;


k = 20;

intervals = [0 0.01 0.1 0.5 1;0 0.01 0.1 0.5 1];
%intervals = [0 0.01 0.1 0.5 1 2 5;0 0.01 0.1 0.5 1 2 5];

noisesigma = [2,4];

potential = @L1;

mu1 = [-1;0]; mu2 = [1,0];
s1 = [0.1,0.1]; s2 = [0.1,0.1];
r1 = mvnrnd(mu1,s1,100); r2 = mvnrnd(mu2,s2,100); r = [r1',r2']';

%noise = mvnrnd([0,0],noisesigma,k); 
lp = laprnd(k,2); lp(:,1) = lp(:,1)*noisesigma(1); lp(:,2) = lp(:,2)*noisesigma(2); noise = lp; rn = [r1',r2',noise']'; 
h = plot(noise(:,1),noise(:,2),'rx','LineWidth',1); hold on; 
plot(r1(:,1),r1(:,2),'ko','LineWidth',1); axis equal; set(gca,'FontSize',15);
plot(r2(:,1),r2(:,2),'ks','LineWidth',1); axis equal; set(gca,'FontSize',15);

figure;

clear pc1L1; 
clear pc1L1opt; 
clear pc1L2; 
clear pc1PQSQL2; 
clear pc1PQSQLSQRT;
clear nps;
ii = 1;
for i=0:step:noisepoints  
display(sprintf('%i',i));     
for j=1:nsamples 
        %noise = mvnrnd([0,0],noisesigma,i); rn = [r1',r2',noise']'; 
        lp = laprnd(i,2); lp(:,1) = lp(:,1)*noisesigma(1); lp(:,2) = lp(:,2)*noisesigma(2); noise = lp; rn = [r1',r2',noise']'; 
        pc = pca(rn); pc1L2(ii,j)=abs(pc(1,1)); 
        pcPQSQ = pcaPQSQ_fast(rn,1,'potential',@L1,'optimize',0,'intervals',intervals); pc1L1(ii,j)=abs(pcPQSQ(1,1)); 
        %pcPQSQopt = pcaPQSQ_fast(rn,1,'potential',@L1,'optimize',2,'intervals',intervals); pc1L1opt(ii,j)=pcPQSQopt(1,1); 
        %pcPQSQL2 = pcaPQSQ_fast(rn,1,'potential',@L2,'optimize',0,'intervals',intervals); pc1PQSQL2(ii,j)=pcPQSQL2(1,1); 
        %pcPQSQLSQRT = pcaPQSQ_fast(rn,1,'potential',@LSQRT,'optimize',0,'intervals',intervals); pc1PQSQLSQRT(ii,j)=pcPQSQLSQRT(1,1); 
    end; 
nps(ii) = i;
ii=ii+1;    
end; 

plot(nps,[mean(pc1L1')],'r-','LineWidth',2); hold on; plot(nps,[mean(pc1L1')-mad(pc1L1')],'r--','LineWidth',1); plot(nps,[mean(pc1L1')+mad(pc1L1')],'r--','LineWidth',1);  
%plot([mean(pc1L1opt')],'m-','LineWidth',2); hold on; plot([mean(pc1L1opt')-mad(pc1L1opt')],'m--','LineWidth',1); plot([mean(pc1L1opt')+mad(pc1L1opt')],'m--','LineWidth',1);  
plot(nps,mean(pc1L2'),'b-','LineWidth',2); plot(nps,mean(pc1L2')-mad(pc1L2'),'b--','LineWidth',1);   plot(nps,mean(pc1L2')+mad(pc1L2'),'b--','LineWidth',1);  
%plot(mean(pc1PQSQL2'),'g-','LineWidth',2); plot(mean(pc1PQSQL2')-mad(pc1PQSQL2'),'g--','LineWidth',1);   plot(mean(pc1PQSQL2')+mad(pc1PQSQL2'),'g--','LineWidth',1);  
%plot([mean(pc1PQSQLSQRT')],'k-','LineWidth',2); hold on; plot([mean(pc1PQSQLSQRT')-mad(pc1PQSQLSQRT')],'k--','LineWidth',1); plot([mean(pc1PQSQLSQRT')+mad(pc1PQSQLSQRT')],'k--','LineWidth',1);  

sumdiff = sum(abs(mean(pc1L1')-mean(pc1L2')))

	