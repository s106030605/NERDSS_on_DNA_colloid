

The Matlab codes which were used to determine the Core-Shell completion from the images is coreshellclusters.m. The code "coreshellanalysis.m", "get\_numerical\_input.m", "get\_textual\_input" are helper code.



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

    % - get\_numerical\_input.m

    % - get\_textual\_input.m

    % - coreshellanalysis.m

