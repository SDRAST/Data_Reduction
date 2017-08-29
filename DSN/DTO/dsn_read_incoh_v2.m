filepath = '/home/kuiper/mnt/dto/home/kuiper/mnt/gpu1/var/tmp/'
#filename = '20160703215524_1'; % almost reasonable quantization
#filename = '20160705174750_1'
#filename = '20160705180028_1'
#filename = '20160705190508_1'
filename = '20160705194100_1'

fd = fopen(strcat(filepath,filename),'r');

end_packet = 2;
data1 = zeros(1,1024);
data2 = zeros(1,1024);

kurtosis_1 = zeros(1,1024);
kurtosis_2 = zeros(1,1024);

figure;
for j=1:end_packet % both polarizations
    discard = fread(fd,1,'uint64');
    discard = fread(fd,1,'uint64'); % fread(fd,1,'ubit64');    

        data_all=fread(fd,2*4*512*2,'uint8'); %8 so that we can separate the bytes!
        data_all_rs=reshape(data_all,2,4,512,2);
        data_all_packed=data_all_rs(2,:,:,:)+data_all_rs(1,:,:,:)*2^8;
        data_all_permute1=permute(data_all_packed,[2,3,4,1]);
        pow_all =data_all_permute1(1:2,:,:);
        kurt_all=data_all_permute1(3:4,:,:);
        
        pow_rs= [pow_all(1,:,:)  pow_all(2,:,:)];
        kurt_rs=[kurt_all(1,:,:) kurt_all(2,:,:)];
        
        pow_perm= permute(pow_rs,[2,3,1]);
        kurt_perm=permute(kurt_rs,[2,3,1]);
        
        subplot(2,1,2)
        plot(abs(kurt_perm)/2^12)
        subplot(2,1,1)
        semilogy(abs(pow_perm))              
  
end;
