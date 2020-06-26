%% Dynamic penalty DS for constrained optimization 
clear all
format long
N_so=50;% size
run=1;

method=2;
method_init='kent';% 'rnd'or 'kent'
global epk fth  max_epoch t1 t2
global lambda
lambda=2/6;
max_epoch=5000;
epk=1;
FR=zeros(1,15);eps=1e-3;
t1=2;t2=5;

    fth==1  % good! t1=2;t2=6; method=2
        fhandle=@fitness2;nonhandle=@st2;
       low_lmt=[10, 40, 90, 270, 270, 270, 330, 330, 355, 0, 0, 0, 600, 6000, 10000];
        up_lmt=[20, 70, 95, 280, 280, 280, 340, 340, 360, 5, 5, 5, 1000, 7000, 10200];
    Dim=length(low_lmt);
        for j=1:run % for every function
            %[x01(j,:),y01(j),y01_con(j,:)]=DSA_mincon_dy(fhandle,nonhandle,method,N_so,Dim,low_lmt,up_lmt,max_epoch,method_init);
             [x01(j,:),y01(j),y01_con(j,:)]=cds_mincon2(fhandle,nonhandle,method,N_so,Dim,low_lmt,up_lmt,max_epoch,method_init);
%         if max(y01_con(j,:)>eps)<1;
%             FR(fth)=FR(fth)+1;
%         end
        end
%         result(fth,:)=[min(y01),mean(y01),max(y01),std(y01)];% for Compare
%         y(fth,:)= y01;  % for calculate Success rateR
     save(['E:\','x01','.dat'], 'x01')
     save(['E:\','y01','.dat'], 'y01')
     save(['E:\','y01_con','.dat'], 'y01_con')



