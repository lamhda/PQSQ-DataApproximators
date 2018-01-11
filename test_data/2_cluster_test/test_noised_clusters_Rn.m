nsamples = 10;
step = 1;
noisepoints = 150;

k = 40;

dimension = 100;

intervals = [0 0.01 0.1 0.5 1];
%intervals = [0 0.01 0.1 0.5 1 2 5];
%intervals = [0 0.01 0.1 0.5 1 2 5;0 0.01 0.1 0.5 1 2 5;0 0.01 0.1 0.5 1 2 5;0 0.01 0.1 0.5 1 2 5];

intervals = repmat(intervals,dimension,1);

potential = @L1;
mu1 = zeros(1,dimension);
mu2 = zeros(1,dimension);

mu1(1) = -1; mu2(1) = 1;

s1 = repmat(0.1,1,dimension); % s1(1) = 0.1; s1(2) = 0.1; s1(5) = 0.1; s1(6) = 0.1;
s2 = repmat(0.1,1,dimension); % s2(1) = 0.1; s2(2) = 0.1; s2(5) = 0.1; s2(6) = 0.1;


r1 = mvnrnd(mu1,s1,100); r2 = mvnrnd(mu2,s2,100); r = [r1',r2']';

noisesigma = repmat(1,1,dimension);
noisesigma(3) = 2;
noisesigma(4) = 4;

noisemu = zeros(1,dimension);
noise = mvnrnd(noisemu,noisesigma,k); 
h = plot(noise(:,1),noise(:,2),'rx','LineWidth',2); hold on; 
plot(r(:,1),r(:,2),'ko','LineWidth',2); axis equal; set(gca,'FontSize',15);

figure;

clear pc1L1; 
clear pc1L1opt; 
clear pc1L2; 
clear pc1PQSQL2; 
clear pc1PQSQLSQRT;
clear u;
clear uPQSQ;
clear nps;
ii = 1;
for i=0:step:noisepoints
display(sprintf('%i',i)); 
for j=1:nsamples 
        %noise = mvnrnd(noisemu,noisesigma,i); rn = [r1',r2',noise']'; 
        lp = laprnd(i,dimension); 
        for dim=1:dimension 
            lp(:,dim) = lp(:,dim)*noisesigma(dim); 
        end
        noise = lp; rn = [r1',r2',noise']';
        
        [pc,u] = pca(rn); pc1L2(ii,j)=computeSimpleTtest(u(:,1:2),100); %pc1L2(ii,j)=pc(1,1); 
        [pcPQSQ,uPQSQ] = pcaPQSQ_fast(rn,2,'potential',@L1,'optimize',0,'intervals',intervals); pc1L1(ii,j)=computeSimpleTtest(uPQSQ(:,1:2),100); %pc1L1(ii,j)=pcPQSQ(1,1); 
        %pcPQSQopt = pcaPQSQ_fast(rn,1,'potential',@L1,'optimize',2,'intervals',intervals); pc1L1opt(ii,j)=pcPQSQopt(1,1); 
        %pcPQSQL2 = pcaPQSQ_fast(rn,1,'potential',@L2,'optimize',0,'intervals',intervals); pc1PQSQL2(ii,j)=pcPQSQL2(1,1); 
        %pcPQSQLSQRT = pcaPQSQ_fast(rn,1,'potential',@LSQRT,'optimize',0,'intervals',intervals); pc1PQSQLSQRT(ii,j)=pcPQSQLSQRT(1,1); 
    end; 
nps(ii) = i;    
ii=ii+1;  
end; 

[v,u] = pca(rn); 
plot(u(201:end,1),u(201:end,2),'rx'); hold on; 
plot(u(1:100,1),u(1:100,2),'ko'); axis equal;
plot(u(101:200,1),u(101:200,2),'ks'); axis equal;
 
figure;

[v,u] = pcaPQSQ_fast(rn,2,'potential',@L1,'intervals',intervals,'optimize',0); 
plot(u(201:end,1),u(201:end,2),'rx','LineWidth',1); hold on; 
plot(u(1:100,1),u(1:100,2),'ko','LineWidth',1); axis equal;
plot(u(101:200,1),u(101:200,2),'ks','LineWidth',1); axis equal;

figure;

plot(nps,[mean(pc1L1')],'r-','LineWidth',2); hold on; plot(nps,[mean(pc1L1')-mad(pc1L1')],'r--','LineWidth',1); plot(nps,[mean(pc1L1')+mad(pc1L1')],'r--','LineWidth',1);  
%plot([mean(pc1L1opt')],'m-','LineWidth',2); hold on; plot([mean(pc1L1opt')-mad(pc1L1opt')],'m--','LineWidth',1); plot([mean(pc1L1opt')+mad(pc1L1opt')],'m--','LineWidth',1);  
plot(nps,mean(pc1L2'),'b-','LineWidth',2); plot(nps,mean(pc1L2')-mad(pc1L2'),'b--','LineWidth',1);   plot(nps,mean(pc1L2')+mad(pc1L2'),'b--','LineWidth',1);  
%plot(mean(pc1PQSQL2'),'g-','LineWidth',2); plot(mean(pc1PQSQL2')-mad(pc1PQSQL2'),'g--','LineWidth',1);   plot(mean(pc1PQSQL2')+mad(pc1PQSQL2'),'g--','LineWidth',1);  
%plot([mean(pc1PQSQLSQRT')],'k-','LineWidth',2); hold on; plot([mean(pc1PQSQLSQRT')-mad(pc1PQSQLSQRT')],'k--','LineWidth',1); plot([mean(pc1PQSQLSQRT')+mad(pc1PQSQLSQRT')],'k--','LineWidth',1);  

sumdiff = sum(abs(mean(pc1L1')-mean(pc1L2')))

%ylim([-0.05, 1.2]);



