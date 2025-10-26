classdef coreshellclusters
    %Object for the analysis of core-shell clusters aggregated in droplets
    %Assumes the input format is 3 tiff files for Channel 1, Channel 2, and Channel 3
    %Channel 1 and 2 are of the fluorescent channels in which particles that form aggregates are shown
    %
    %This code calculates the degree of core-shellness of multi-component,
    %spherical objects.
    %The two components need to be measured independently in two channels
    %of the microscopy image. This code takes two images of the same object
    %in different channels and compares them. The outputs are
    % - an azimuthal average of the intensity in channel 1
    % - an azimuthal average of the intensity of channel 2
    % - the Jensen-Shannon divergence: a measure of how different the two
    % distributions (= azimuthal averages) are (i.e. the number of bits required to transform one
    % distribution to the other, compared to the number of bits required to
    % create the average distribution)
    % - the distance between the median of the two distributions.
    %By convenstion, channel 1 is blue and channel 2 is red
    %
    %Requires helper codes:
    % - get_numerical_input.m
    % - get_textual_input.m
    % - coreshellanalysis.m
    %
    %Shorthand user manual
    %-Choose a name for your experiment, here obj
    %-Initiate an object of type coreshelclusters by typing
    %obj=coreshellclusters()
    %You will be prompted to select an image stack of type tiff.
    %This stack should be 2 images and contain the two channels of your
    %experiment
    %-Analyze yuor data by running obj=analyze(obj,val1,val2,val3)
    %Val1,2,3 need to be numerical values that are chosen appropriately,
    %otherwise the code will not work as intended.
    %- val1 is a blur_length and indicates over how many pixels the image
    %gets blurred. This needs to be chosen such that the clusters in the
    %image become one continuous shape after blurring, but so that two
    %separate clusters stay two discontinuous shapes. So choose the blur
    %length smaller than half the cluster size in pixels but bigger than
    %the size of individual particles in the cluster
    %- val2 is the size of an image that contains one entire cluster, but
    %no more than one cluster, so that the image can be segmented into
    %subimages with each one cluster in it. Choose val 2 larger than the
    %clustersize in pixels, but smaller than the distance between pixels.
    %- val2 is a threshold, H, in the intensity of the image. Any value
    %below H is considered background. Choose H higher than the background
    %intensity, but much smaller than the intensity of the clusters.
    %
    %- Your experiment is analyzed and you can now plot the values of
    %interest using
    % - plot_radial_average(obj)
    % - obj=calc_av_JS_divergence(obj)
    % - obj=calc_av_dw(obj)
    % All calculated values are saved in the object with the name you chose
    %(here obj) in the workspace. Note that each time you run a code, this
    %object overrides. To avoid overriding, just chose a new name:
    %obj2=calc_av_dw(obj) will copy obj and call it obj2 and calculate dw, the
    %average distance between the medians of the distributions, in obj2
    %Have fun!
    %
    %Modification history
    %July 2023: created by Pepijn Moerman
    %
    properties
        fname       %List of 2 filenames, corresponding to 2 channels
        ch1         %The image of channel 1
        ch2         %The image of channel 2
        ch1max
        ch2max
        combined    %The two channels added together
        size        %image size in [width pixels height pixels]
        pos         %list of center of mass of clusters
        improcparams%List of image processing parameters
        cc          %Cluster object
        blocksize   %edgelength of image that contains one cluster
        rejectedlist %list of rejected objects and reason
        clusterlist %list of objects of type coreshellanalysis
        JS_list     %list of JS divergence values for each cluster
        dw_list     %list of distances between red and blue cluster edge.
        
    end
    methods
        function obj=coreshellclusters()
            %connect the object to a tiffstack. Save file and foldername
            %!!! IT IS ESSENTIAL THAT THE THREE FILES HAVE THE SAME NAME
            %AND ONLY THE LAST NUMBER DIFFERS, BEING 1 FOR CHANNEL 1, 2 FOR
            %2 ETC. !!!
            disp('Select the folder that a tiff of channel 1')
            [filen, foldername_temp] = uigetfile('*.tif','Select tif stack');
            fname1=[foldername_temp filen];
            
            fname2=[foldername_temp filen(1:end-5) '2.tif'];
            try A = imread(fname2);
            catch
                disp(['The file with name ' fname2 ' could not be opened'])
                error('Make sure your files are labeled "filename1.tif", "filename2.tif", etc for the 3 channels')
            end
            
            obj.fname={fname1 fname2};
            obj.ch1=imread(fname1);
            obj.ch2=imread(fname2);
            obj.combined=obj.ch1+obj.ch2;
            obj.size=size(obj.ch1);
        end
        function obj=analyze(obj,blur_length,blocksize,H)
            if nargin<4
                H=100;
            end
            if nargin<3
                blocksize=700;
            end
            if nargin<2
                blur_length=20;
            end
            obj=find_clusters(obj,blur_length);
            obj=export_clusters(obj,blocksize);
            obj=analyze_clusters(obj,H);
        end
        function plot_image(obj,i)
            %i is channel number to plot
            %i=3 means combined image
            if nargin<2
                i=1; %default is image 1
            end
            if i>3
                error('you only have two channels and a joint one. Choose i=[1,2,3]')
            end
            if i==1
                A=obj.ch1;
            elseif i==2
                A=obj.ch2;
            elseif i==3
                A=obj.combined;
            end
            figure
            hold on
            colormap('gray'),imagesc(A);
            axis([0 obj.size(2) 0 obj.size(1)])
            hold off
        end
        function plot_image_max(obj,i)
            %i is channel number to plot
            %i=3 means combined image
            if nargin<2
                i=1; %default is image 1
            end
            if i>2
                error('you only have two channels and a joint one. Choose i=[1,2,3]')
            end
            if i==1
                A=obj.ch1max;
            elseif i==2
                A=obj.ch2max;
            end
            figure
            hold on
            colormap('gray'),imagesc(A);
            axis([0 obj.size(2) 0 obj.size(1)])
            hold off
        end
        function obj=maximize_contrast(obj)
            im1=double(obj.ch1);
            maxval1=max(max(im1));
            minval1=min(min(im1));
            im1=(im1-minval1)./(maxval1-minval1);
            obj.ch1max=uint16(im1.*(2^16-1));
            
            im2=double(obj.ch2);
            maxval2=max(max(im2));
            minval2=min(min(im2));
            im2=(im2-minval2)./(maxval2-minval2);
            obj.ch2max=uint16(im2.*(2^16-1));
        end
        function plot_bothchannels(obj)
            im1=obj.ch1max;
            im2=obj.ch2max;
            figure
            %define colormaps
            redcolormap=[linspace(0,1,256)', zeros(256,2)];
            cyancolormap=[zeros(256,1),linspace(0,1,256)',linspace(0,1,256)'];
            %plot first data 
            ax1 = axes;
            im = imagesc(ax1,im1);
            im.AlphaData = 0.5; % change this value to change the background image transparency
            hold all;
            %plot second data
            ax2 = axes;
            im1 = imagesc(ax2,im2);
            im1.AlphaData = 0.5; % change this value to change the foreground image transparency
            %link axes
            linkaxes([ax1,ax2])
            %%Hide the top axes
            ax2.Visible = 'off';
            ax2.XTick = [];
            ax2.YTick = [];
            %add differenct colormap to different data if you wish
            colormap(ax1,cyancolormap)
            colormap(ax2,redcolormap)
        end
        function obj=find_clusters(obj,blur_length,threshold)
            if nargin<2
                blur_length=2;
            end
            if nargin<3
                threshold=0.05;
            end
            check='n';
            while check~='y'
                figure
                A=imgaussfilt(obj.combined,blur_length);
                colormap('gray'),imagesc(A);
                axis([0 obj.size(2) 0 obj.size(1)]);
                message=('Is this a good blurring lengthscale(y/n): ');
                acceptables=['y','n'];
                check=get_textual_input(message, acceptables);
                if check=='n'
                    blur_length=get_numerical_input('What length scale do you want to use instead: ');
                end
            end
            
            meanA=mean(mean(A));
            A(find(A<threshold*meanA))=0;
           
            obj.improcparams={'Blur length' blur_length ; 'Threshold' threshold};
            CC = bwconncomp(A);
            obj.cc=CC; 
            obj.pos=NaN(CC.NumObjects,2);
            
            for i=1:CC.NumObjects
                Ai=zeros(obj.size(1),obj.size(2));
                Ai(CC.PixelIdxList{1,i})=1;
                [c, r] = find(Ai == 1);
                obj.pos(i,:) = [mean(r), mean(c)];
            end
            close all
            
            figure
            hold on
            box on
            colormap('gray'),imagesc(obj.combined);
            viscircles(obj.pos(:,1:2),ones(CC.NumObjects,1)*blur_length);
            axis([0 obj.size(2) 0 obj.size(1)]);
            hold off

        end
        function obj=export_clusters(obj,blocksize)
            if isempty(obj.ch1max)
                obj=maximize_contrast(obj);
            end
            obj.clusterlist={};
            obj.rejectedlist={};
            if nargin<2
                if isempty(obj.blocksize)
                    obj.blocksize=get_numerical_input('What is the size of the clusters in pixels: ');
                end
            else
                obj.blocksize=blocksize;
            end
%             if nargin<3
%                 H=100;
%             end
            
            edge=round(obj.blocksize/2);
            count=0;
            count2=0;
            for i=1:length(obj.pos(:,1))
                xi=round(obj.pos(i,2));
                yi=round(obj.pos(i,1));
                if xi-edge>0 && xi+edge<obj.size(1) && yi-edge>0 && yi+edge<obj.size(2)
                    A_ch1=obj.ch1max(xi-edge:xi+edge,yi-edge:yi+edge);
                    A_ch2=obj.ch2max(xi-edge:xi+edge,yi-edge:yi+edge);
                    ai=coreshellanalysis(A_ch1,A_ch2);
                    ai=azimuthalaverage(ai);
                    plot_azimuthal_average(ai)
                    message='Is this cluster acceptable (y/n): ';
                    acceptables=['y', 'n'];
                    answer=get_textual_input(message,acceptables);
                    options={'Size out of range', 'Intensity too low', 'Other'};
                    if answer=='y'
                        count=count+1;
                        obj.clusterlist{count}=ai;
                    else
                        moveon=0;
                        while moveon==0
                            disp('Reason for rejection (1/2/3): ')
                            message2='1) Size out of range. 2) Intensity too low. 3) Other.';
                            answer2=get_numerical_input(message2);
                            if answer2==1 || answer2==2 || answer2==3
                                count2=count2+1;
                                obj.rejectedlist{count2,1}=ai;
                                obj.rejectedlist{count2,2}=options{answer2};
                                moveon=1;
                            end
                        end
                    end
                    
                end
            end
            close all;
        end
        function obj=analyze_clusters(obj,H)
            if nargin<2
                H=100;
            end
            for i=1:length(obj.clusterlist)
                ai=obj.clusterlist{i};
                ai=calc_JS_divergence(ai);
                ai=calc_shell_width(ai,H);
                obj.clusterlist{i}=ai;
            end
        end
        function plot_radial_averages(obj)
            try ai=obj.clusterlist{1};
            catch 
                error('You have not identified any cluters. Run export_clusters first')
            end
            if isempty(ai.H)
                error('You have not analyzed the clusters. Run analyze_clusters first')
            end
            %channel 1 is always blue.
            %channel 2 is always red
            figure
            hold on
            box on
            y_ch1=[];
            y_ch2=[];
            xtemp=[];
            xall=linspace(0,1.2,200)';
            for i=1:length(obj.clusterlist)
                ai=obj.clusterlist{1,i};
                x = ai.azav_ch1(:,1);
                xnonzero_ch1=find(ai.azav_ch1(:,2)>ai.H);
                xnow=x./xnonzero_ch1(end);
                %xtemp=[xtemp, xnow];
                ych1new(:,1)=interp1(xnow', ai.azav_ch1(:,2)', xall', 'linear','extrap')';
                ych2new(:,1)=interp1(xnow', ai.azav_ch2(:,2)', xall', 'linear','extrap')';
                ych1sum=sum(ych1new(:,1));
                ych2sum=sum(ych2new(:,1));
                y_ch1=[y_ch1 , ych1new(:,1)./ych1sum];
                y_ch2=[y_ch2 , ych2new(:,1)./ych2sum];
                plot(xnow, ai.azav_ch1(:,2),'b.')
                plot(xnow, ai.azav_ch2(:,2),'r.')
                %plot(xall, ych1new(:,1),'b.')
                %plot(xall, ych2new(:,1),'r.')
                xlabel('r (pixels)')
                ylabel('Intensity (counts)')
            end
            hold off
            
            
            ymeanch1=mean(y_ch1,2);
            ymeanch2=mean(y_ch2,2);
            
            figure
            hold on
            box on
            plot(xall, ymeanch1(:,1),'b.')
            plot(xall, ymeanch2(:,1),'r.')
            
            curve1 = mean(y_ch1,2)+std(y_ch1',1)';
            curve2 = mean(y_ch1,2)-std(y_ch1',1)';
            plot(xall, curve1, 'b-', 'LineWidth', 1);
            plot(xall, curve2, 'b-', 'LineWidth', 1);
            patch([xall' fliplr(xall')], [curve1' fliplr(curve2')], 'b','FaceAlpha',0.3)
            
            curve3 = mean(y_ch2,2)+std(y_ch2',1)';
            curve4 = mean(y_ch2,2)-std(y_ch2',1)';
            plot(xall, curve3, 'r-', 'LineWidth', 1);
            plot(xall, curve4, 'r-', 'LineWidth', 1);
            patch([xall' fliplr(xall')], [curve3' fliplr(curve4')], 'r','FaceAlpha',0.3)
            pbaspect([1 1 1])
            %axis([0,1.2,0,max([max(curve3) max(curve1)])])
            %axis([0,1.2,0,1.4*max([max(ymeanch1) max(ymeanch2)])])
            axis([0,1.2,0,0.015])
            hold off
        end
        function plot_radial_average_i(obj,i)
            try ai=obj.clusterlist{1};
            catch 
                error('You have not identified any cluters. Run export_clusters first')
            end
            if isempty(ai.H)
                error('You have not analyzed the clusters. Run analyze_clusters first')
            end
             ai=obj.clusterlist{1,i};
             plot_azimuthal_average(ai);
        end
        function obj=calc_av_JS_divergence(obj)
            try ai=obj.clusterlist{1};
            catch 
                error('You have not identified any cluters. Run export_clusters first')
            end
            if isempty(ai.H)
                error('You have not analyzed the clusters. Run analyze_clusters first')
            end
                %output= [av std]
                JSlisttemp=NaN(length(obj.clusterlist),1);
                for i=1:length(obj.clusterlist)
                	a=obj.clusterlist{i};
                    JSlisttemp(i,1)=a.JS;
                end
                obj.JS_list=JSlisttemp;
                
                figure
                hold on
                box on
                histogram(obj.JS_list)
                hold off
                
                disp([mean(JSlisttemp) std(JSlisttemp)]);
        end
        function obj=calc_av_dw(obj)
            try ai=obj.clusterlist{1};
            catch 
                error('You have not identified any cluters. Run export_clusters first')
            end
            if isempty(ai.H)
                error('You have not analyzed the clusters. Run analyze_clusters first')
            end
            dwlisttemp=NaN(length(obj.clusterlist),1);
            for i=1:length(obj.clusterlist)
                a=obj.clusterlist{i};
                dwlisttemp(i,1)=a.dw_shell;
            end
            obj.dw_list=dwlisttemp;
            
            figure
            hold on
            box on
            histogram(obj.dw_list)
            hold off
            
            disp([mean(dwlisttemp,'omitnan') std(dwlisttemp,'omitnan')]);
        end
    end%end methods
end%end class