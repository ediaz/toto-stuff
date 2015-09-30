		fide  = fopen('/var/tmp/green1.rsf@','r'); 
		temp = fread(fide,[166,101*476],'float32');
		fclose(fide);
		%gp(:,ix,:) = temp(1:ntf,:);
        Green=reshape(temp,166,101,476);
        close all
        
        for i=1:10:300
            subplot(2,1,1)
            imagesc(Green(1:91,:,i));
            subplot(2,1,2)
            imagesc(squeeze(gfull(i*2,:,:,3))');
            pause(0.5);
         end     