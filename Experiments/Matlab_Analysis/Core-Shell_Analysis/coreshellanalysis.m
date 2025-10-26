classdef coreshellanalysis
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
    %Operation manual:
    % 
    %
    %Modification history:
    %   2022: Created by Pepijn Moerman

    %properties
    %    ch1         %the cluster image in channel 1
    %    ch2         %the cluster image in channel 2
    %    rtheta      %r and theta values of grid
    %    size        %the cluster image size
    %    azav_ch1    %azimuthal average of channel 1
    %    azav_ch2    %azimuthal average of channel 2
    %    JS          %Jenson-Shannon divergence of the two color distributions
    %    H           %height in counts, used to find dw_shell in pixels
    %    dw_shell    %widht of the shell in pixels (CHANNEL 1 - CHANNEL 2)
    %end

    properties
        ch1         % Image in channel 1 (e.g., blue fluorescence)
        ch2         % Image in channel 2 (e.g., red fluorescence)
        rtheta      % Radial and angular coordinates of the grid
        size        % Image size (rows, columns)
        azav_ch1    % Azimuthal average of intensity in channel 1
        azav_ch2    % Azimuthal average of intensity in channel 2
        JS          % Jensen-Shannon divergence (difference between distributions)
        H           % Height threshold for shell width calculation
        dw_shell    % Difference in shell width (pixels) between channel 1 and 2
    end

    methods
        function obj=coreshellanalysis(ch1,ch2)
            obj.ch1=ch1;
            obj.ch2=ch2;
            obj.size=size(ch1);
            
            %calculate the r and theta values of each point in the grid
            obj.rtheta=NaN(obj.size(1),obj.size(2),2);
            cx=floor(obj.size(2)/2)+1;
            cy=floor(obj.size(1)/2)+1;
            for x=1:obj.size(2)
                for y=1:obj.size(1)
                    obj.rtheta(y,x,1)=sqrt((x-cx)^2+(y-cy)^2); %Calculate the radius
                    obj.rtheta(y,x,2)=atan2(y-cy,x-cx);
                end
            end
        end
        function obj=azimuthalaverage(obj,dr)
            if nargin<2
                dr=1;
            end
            dr=round(dr);
            rmax=round(obj.size(1)/2);
            rvals=obj.rtheta(:,:,1);
            count=0;
            nmax=round((rmax-dr)/dr);
            av_ch1=NaN(nmax,2);
            av_ch2=NaN(nmax,1);
            %loop over all r
            for r=1:dr:rmax-dr
                count=count+1;
                rlow=r;
                rhigh=r+dr;
                indices=find(rvals(:,:)>rlow & rvals(:,:)<rhigh);
                av_ch1(count,1)=r;
                av_ch1(count,2)=mean(obj.ch1(indices));
                av_ch2(count,1)=r;
                av_ch2(count,2)=mean(obj.ch2(indices));
            end
            obj.azav_ch1=av_ch1-min(av_ch1);
            obj.azav_ch2=av_ch2-min(av_ch2);
            
        end
        function plot_azimuthal_average(obj)
            %channel 1 is always blue.
            %channel 2 is always red
            figure
            hold on
            box on
            plot(obj.azav_ch1(:,1), obj.azav_ch1(:,2),'b.')
            plot(obj.azav_ch2(:,1), obj.azav_ch2(:,2),'r.')
            xlabel('r (pixels)')
            ylabel('Intensity (counts)')
            hold off
        end
        function obj=calc_JS_divergence(obj)
            %calculates and outputs the JS divergence of the two color
            %distributions in this object (azav_ch1 and azav_ch2)
            
            %normalize azav_ch1 and azav_ch2 to themselves
            dch1=obj.azav_ch1(:,2)./sum(obj.azav_ch1(:,2));
            dch2=obj.azav_ch2(:,2)./sum(obj.azav_ch2(:,2));
            
            %calculate the mean of the two distributions
            M=mean([dch1(:,1) dch2(:,1)],2);
            
            %calculate the Kullbeck-Leibler divergence for both
            %distributions
            KLprep_ch1=NaN(length(M(:,1)),1);
            KLprep_ch2=NaN(length(M(:,1)),1);
            for i=1:length(M(:,1))
                KLprep_ch1(i,1)=dch1(i,1)*log2(dch1(i,1)/M(i,1));
                KLprep_ch2(i,1)=dch2(i,1)*log2(dch2(i,1)/M(i,1));
            end
            KL_ch1=sum(KLprep_ch1,'omitnan');
            KL_ch2=sum(KLprep_ch2,'omitnan');
            
            %Calculate the JS divergence
            obj.JS=(KL_ch1+KL_ch2)/2;
        end
        function obj=calc_shell_width(obj,H)
            %Define the length scale of particle type as the median of the
            %density distribution
            if nargin<2
                H=100;
            end
            
            %find the integral of both distributions
            ych1sum=sum(obj.azav_ch1(:,2));
            ych2sum=sum(obj.azav_ch2(:,2));
            
            index1=find(cumsum(obj.azav_ch1(:,2))>0.5*ych1sum)
            index2=find(cumsum(obj.azav_ch2(:,2))>0.5*ych2sum)
            obj.dw_shell=index1(1,1)-index2(1,1);
            
%             figure
%             hold on
%             box on
%             plot(obj.azav_ch1(:,1), obj.azav_ch1(:,2),'b.')
%             plot([index1(1,1) index1(2,1)]',[0 10^4]','b-')
%             plot(obj.azav_ch2(:,1), obj.azav_ch2(:,2),'r.')
%             plot([index2(1,1) index2(2,1)]',[0 10^4]','r')
%             xlabel('r (pixels)')
%             ylabel('Intensity (counts)')
%             hold off
            
%             dwtemp=NaN(H,1);
%             for h=1:H
%                 %find the values above height H for both channels
%                 indicesch1=find(obj.azav_ch1(:,2)>h);
%                 indicesch2=find(obj.azav_ch2(:,2)>h);
%                 
%                 %subtract the position to find dw
%                 if ~isempty(indicesch1) && ~isempty(indicesch2)
%                     endval1=length(indicesch1);
%                     endval2=length(indicesch2);
%                     dwtemp(h,1)=indicesch1(endval1)-indicesch2(endval2);
%                 else
%                     dwtemp(h,1)=NaN;
%                 end
%             end
            %save values
            obj.H=H;
%             obj.dw_shell=mean(dwtemp(:,1),'omitnan');
        end
    end
end
