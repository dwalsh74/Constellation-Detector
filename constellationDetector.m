
clear all;
close all;
clc;

%% Parameters

%Note: Make sure the input image and template images are included in the
%      current MATLAB path

%Enter the test image filename below
inputImStr = ('Gemini_Googled.png');

TemplateList= {'Aries','Taurus','Gemini','Cancer','Leo','Libra',...
                'Capricorn','Aquarius','Pisces','Virgo'};

%Configure Figure Font
set(0,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%Dictates whether noise and random rotation will be applied to the image
distort = false;

%How much the size can deviate and be considered a match
scaleTolerance = 0.18;

%How much the distance can deviate and be considered a match
distTolerance = 0.18;

%How many degrees the angle can deviate and be considered a match
angleTolerance = 15; 

%Run loop 'runs' number of times 
runs = 1;

%% Main Code
for runIndex = 1:runs

    x=length(TemplateList); 
    for tempnum = 1:1:length(TemplateList)
        
        %% Template Conversion
        % Convert template image into numerical data  
        
        %Get the current template name
        TemplateName = TemplateList{tempnum};
        
        %Output Progress to console
        msg = sprintf('Searching for %s constellation...',TemplateName);
        disp(msg);

        %Create Image names based on DesiredConstellation Global
        templateStr = append(TemplateName,'_Template.png');
        templateLinesStr = append(TemplateName,'_Template_Lines.png');


        %Read template images
        binaryTemplate = imread(templateStr);

        %Convert to grayscale if is color image
        if size(binaryTemplate,3) == 3
            binaryTemplate = rgb2gray(binaryTemplate);
        end
        binaryTemplate_Lines = imread(templateLinesStr) > 0;
        
        %Get image dimensions for padding
        heightPad = length(binaryTemplate(:,1))*2;
        widthPad = length(binaryTemplate(1,:))*2;
        
        %Add padding to accomodate rotation later in the code
        binaryTemplate = padarray(binaryTemplate,...
                                    [widthPad heightPad],0,'both');
        binaryTemplate_Lines = padarray(binaryTemplate_Lines,...
                                        [widthPad heightPad],0,'both');

        % Binarize Image
        binaryTemplate = binaryTemplate > 0;

        binaryTemplate_Lines =  binaryTemplate_Lines - binaryTemplate;

        %Get stars as objects
        [L,num] = bwlabel(binaryTemplate);

        %Get Template Centers and Areas
        topn = bwareafilt(binaryTemplate, num, 'Largest');

        %Get area and position of every star in the template
        Regions = regionprops(topn,'Area','Centroid');
        templateCenters_Temp = [Regions.Centroid];
        templateAreas = [Regions.Area];

        templateCenters = zeros(num,2);
        j = 1;
        for i = 1:2:num*2
            templateCenters(j,1) = templateCenters_Temp(i);
            templateCenters(j,2) = templateCenters_Temp(i+1);
            j = j+1;
        end

        %Use bubble sort to rank stars in template by size
        while true
            wasChanged = 0;

            for i = 1:(length(templateAreas)-1)

                % Last i elements are already in place
                if (templateAreas(i) < templateAreas(i+1))

                    %Swap values in templateAreas
                    temp1 = templateAreas(i);
                    temp2 = templateAreas(i+1);
                    templateAreas(i) = temp2;
                    templateAreas(i+1) = temp1;

                    %Swap values in templateCenters
                    temp1 = templateCenters(i,1);
                    temp2 = templateCenters(i,2);
                    temp3 = templateCenters(i+1,1);
                    temp4 = templateCenters(i+1,2);
                    templateCenters(i,1) = temp3;
                    templateCenters(i,2) = temp4;   
                    templateCenters(i+1,1) = temp1;
                    templateCenters(i+1,2) = temp2;                

                    %Indicate a value was changed
                    wasChanged = 1;
                end
            end

            % Exit loop if nothign was changed this iteration
            if wasChanged == 0
               break
            end

        end    

        %Get parameters of the largest star in the template image
        largestTemplateStar = bwareafilt(binaryTemplate, 1, 'Largest');
        Regions = regionprops(largestTemplateStar,'Area','Centroid');
        largestTempArea = Regions.Area;
        largestTempCent = [Regions.Centroid];

        %Get distances and size ratios of each star wrt the largest areas
        templateScaleFactors = zeros(num,1);
        templateDists = zeros(num,2);

        for i = 1:num

           templateScaleFactors(i) = templateAreas(i) / largestTempArea;
           templateDists(i,1) = largestTempCent(1) - templateCenters(i,1);
           templateDists(i,2) = largestTempCent(2) - templateCenters(i,2);

        end

        %% Input Image Processing
        %  Process the input image and get numerical data from stars
        
        %Read input image
        if(tempnum == 1)
            inputIm = imread(inputImStr);

            %Add noise and random rotation if enabled in parameters
            if (distort == 1)
                %Rotate input image by a random angle
                rotationConstant = randi(360);
                inputIm = imrotate(inputIm,rotationConstant);

                %Add noise to input image
                inputIm = imnoise(inputIm,'salt & pepper',0.01);
            end

            %Save original image to display at end
            originalIm = inputIm;

            %Add whitespace to right side of originalIm
            originalIm = padarray(originalIm,[0 50],255,'both');

        else       
            inputIm = nextiterationIm;  
        end
        
        im2Gray = rgb2gray(inputIm);
        im2BW = im2Gray > 30 ;

        %Remove extremely small stars
        struct = strel('disk',1);
        im2BW = imopen(im2BW,struct);
        
        %Fill holes
        struct = strel('disk',2);%2
        im2BW = imclose(im2BW,struct);

        %Remove small stars
        struct = strel('square',2);%2
        im2BW = imopen(im2BW,struct);

        templateAngle = 0;
        inputImAngle = 0;

        %Get Centers and areas of largest n stars in the image
        numStars = 500;

        %Don't test more stars than are in the image
        [L1, num1] = bwlabel(im2BW);
        if numStars > num1
            numStars = num1;
        end

        tempIm = im2BW;
        top = im2BW*0;

        warning('off');
        for i = 1:numStars
            currentLargest = bwareafilt(tempIm, 1, 'Largest');
            tempIm = (tempIm - currentLargest) > 0;
            top = (top + currentLargest) > 0;
        end    
        warning('on');

        Regions=regionprops(top,'Area','Centroid');
        actualAreas = [Regions.Area];
        actualCenters_Temp = [Regions.Centroid];
        actualCenters = zeros(numStars,2);

        j = 1;
        for i = 1:2:numStars*2
            actualCenters(j,1) = actualCenters_Temp(i);
            actualCenters(j,2) = actualCenters_Temp(i+1);
            j = j+1;
        end

        centerComp(1,1) = largestTempCent(1);
        centerComp(1,2) = largestTempCent(2);

        %% Star Detection and Scoring
        %Attempt to find the current constellation in the input image
        
        %The max number of stars the algoritm will test as reference stars
        maxTries = 50;

        %Don't test more stars than are in the image
        if maxTries > num1
            maxTries = num1;
        end
        
        %Iterate through potential reference stars, giving each a score
        tempIm = im2BW;
        for n = 1:maxTries

            if n ==2
                afg=1;
            end


            angle = 0;

            actualMinY = 99999;
            actualMaxY = 0;
            actualMinX = 99999;
            actualMaxX = 0;

            numTrue = 0;
            count = 0;

            % Get nth biggest star
            warning('off');
            nthStar = bwareafilt(tempIm, 1, 'Largest');
            tempIm = (tempIm - nthStar) > 0;
            warning('on');

            %Get star properties
            Regions=regionprops(nthStar,'Area','Centroid');
            nthStarArea = Regions.Area;
            nthStarCentroid = [Regions.Centroid];

            scaleFactors = zeros(1,num);
            angles = zeros(1,num);
            
            %Iterate through all stars in the template, attempting to find
            %a match in the input image
            for i = 1:num

                %Get star size ratio
                templateSCF = templateAreas(i)/max(templateAreas);

                %Compare current star center to largest star center
                centerComp(2,1) = templateCenters(i,1);
                centerComp(2,2) = templateCenters(i,2);

                %Skip loop where star is compared against itself
                if centerComp(1,1) == centerComp(2,1)
                    if centerComp(1,2) == centerComp(2,2)
                        continue;
                    end
                end

                templateDist = pdist(centerComp,'euclidean');

                %Get Angle Between Stars in Template
                template_xdiff = centerComp(1,1) - centerComp(2,1);
                template_ydiff = centerComp(1,2) - centerComp(2,2);
                templateAngled = atan2d(template_ydiff,template_xdiff);

                %Correct for reference image rotation
                templateAngled = templateAngled - angle;

                % Correct angles above 180 deg or below -180 deg 
                if templateAngled < -180
                    templateAngled = templateAngled + 360;
                elseif templateAngled > 180
                    templateAngled = templateAngled - 360;
                end

                %Pass in 1st star & size/dist info
                retArray = starCheck(actualCenters,actualAreas,...
                            largestTempArea,templateSCF,...
                            scaleTolerance,distTolerance,templateDist,...
                            nthStarArea,nthStarCentroid,templateAngled,...
                            numTrue,angleTolerance);

                if retArray(1) ~= 0

                    %Template Angles
                    template_xdiff = centerComp(1,1) - centerComp(2,1);
                    template_ydiff = centerComp(1,2) - centerComp(2,2);
                    templateAngled = atan2d(template_ydiff,template_xdiff);

                    nthStarDist = pdist(retArray,'euclidean'); 

                    numTrue = numTrue + 1; 

                    %Input Image Angles
                    xdiff2 = retArray(1,1) - retArray(2,1);
                    ydiff2 = retArray(1,2) - retArray(2,2);
                    inputImAngled = atan2d(ydiff2,xdiff2);

                    if numTrue == 1
                        templateAngle = templateAngled;
                        inputImAngle = inputImAngled;
                        angle = templateAngle - inputImAngle;
                    end

                    %Record scale factor from this pair of stars 
                    scaleFactors(i) = nthStarDist/templateDist;
                    angles(i) = templateAngled - inputImAngled;

                end

            end

            %Exit for loop if match is found
            if num > 5
                passingScore = floor(num*0.7);
            else
                passingScore = floor((num*0.7));
            end

            if numTrue == 12
            xsdf = 1;
            end
            if numTrue > passingScore

                %Clean zeros from scale and rotation factors
                scaleFactors = scaleFactors(scaleFactors ~= 0);
                angles = angles(angles ~= 0);

                %Calculate scale factor for template overlay
                scaleComp = median(scaleFactors);
                angleComp = median(angles);
                break;

            end

        end

        %% Constellation Mapping and Overlay
        %  Format and overlay the constellation onto the image
        
        if numTrue > passingScore
            
            msg = sprintf('%s Constellation Detected!',TemplateName);
            disp(msg);

            %Get angle
            angle = angleComp;

            %Rotate Image
            rotTemplate = imrotate(binaryTemplate,angle);

            %Get X and Y extrema of stars in scaled, shifted template
            Regions = regionprops(rotTemplate,'Centroid');
            rotTemplateCenters_temp = [Regions.Centroid];    

            rotatedTemplateCenters = zeros(num,2);
            j = 1;
            for i = 1:2:num*2
                rotatedTemplateCenters(j,1) = rotTemplateCenters_temp(i);
                rotatedTemplateCenters(j,2) = rotTemplateCenters_temp(i+1);
                j = j+1;
            end

            templateMinY = min(rotatedTemplateCenters(:,2)); 
            templateMaxY = max(rotatedTemplateCenters(:,2)); 
            templateMinX = min(rotatedTemplateCenters(:,1)); 
            templateMaxX = max(rotatedTemplateCenters(:,1)); 

            %Resize template image to match output image based on extrema
            templateDeltaY = templateMaxY - templateMinY;
            templateDeltaX = templateMaxX - templateMinX;

            actualDeltaY = actualMaxY - actualMinY;
            actualDeltaX = actualMaxX - actualMinX;

            scaledTemplate = imresize(rotTemplate,scaleComp);

            top = bwareafilt(scaledTemplate, 1, 'Largest');
            Regions = regionprops(top,'Centroid');
            scaledTemplateCentroid = [Regions.Centroid];

            %Shift template to match coordinates of original
            xdiff = nthStarCentroid(1) - scaledTemplateCentroid(1);
            ydiff = nthStarCentroid(2) - scaledTemplateCentroid(2);

            %Prepare template for overlay over output image
            rotTemplate_Lines = imrotate(binaryTemplate_Lines,angle);
            scaledTemplate_Lines = imresize(rotTemplate_Lines,scaleComp);
            shiftedTemplate_Lines = imtranslate(scaledTemplate_Lines,...
                                    [xdiff, ydiff]) > 0;

            %Overlay (addition for now)
            [x1, y1] = size(shiftedTemplate_Lines);
            [x2, y2] = size(im2BW);

            minWidth = min(x1,x2);
            minHeight = min(y1,y2);

            for i = 1: minWidth
                for j = 1:minHeight
                    if shiftedTemplate_Lines(i,j) > 0
                        inputIm(i,j,1) = 0;
                        inputIm(i,j,2) = 255;
                        inputIm(i,j,3) = 0;
                    end
                end
            end

            %Find largest star in template
            largestStar = bwareafilt(shiftedTemplate_Lines, 1, 'Largest');
            Regions = regionprops(largestStar,'Area','Centroid');
            center = [Regions.Centroid];

            %Add label to image
            try
                textPosition(1) = center(1) + 20;
                textPosition(2) = center(2) + 30;
            catch
                textPosition(1) = center(1) - 20;
                textPosition(2) = center(2) - 30;
            end
            inputIm = insertText(inputIm,textPosition,TemplateName,...
                    'BoxOpacity',0,'TextColor','White','FontSize',36);
        end

        nextiterationIm = inputIm; 

    end

    %% Final Output

    figure;
    imshowpair(originalIm,inputIm,'montage');
    title('Original Input Image Vs. Final Output Image');

end

%% Functions

%Search the input image for a match to a given pair of stars
function retArray = starCheck(inputImCenters,inputImageAreas,...
                                largestTemplateStarArea,tempSCF,...
                                scaleTolerance,distTolerance,...
                                tempDist,nthStarArea,...
                                nthStarCentroid,templateAngled,numTrue,...
                                angleTolerance)
    retArray = zeros(4,1);
    numLoops = length(inputImageAreas);

    bestDistMatch = 0;
    bestScaleMatch = 0;
    bestAngleMatch = 0;

    %Find pairs
    for i = 1:numLoops      

        %Get star center and radii
        currentStarArea = inputImageAreas(i);
        currentStarCentroid(1) = inputImCenters(i,1);
        currentStarCentroid(2) = inputImCenters(i,2);     

        %Get scale factor
        currSCF = currentStarArea/nthStarArea;

        % Check if size ratio of current star pair is within the tolerance
        if (currSCF >= (tempSCF-(tempSCF*scaleTolerance)))
        if (currSCF <= tempSCF+(tempSCF*scaleTolerance))

            %Get distance
            centerCompare = [nthStarCentroid ; ...
                                inputImCenters(i,1) inputImCenters(i,2)];
            currentPairDist = pdist(centerCompare,'euclidean'); 

            %Get scale factor of biggest in template to biggest in pair
            SCFComp = sqrt(nthStarArea) / sqrt(largestTemplateStarArea);

            %Scale distance based on scale factor
            ScaledDist = currentPairDist/SCFComp;

            %Compare distances
            if (ScaledDist >= (tempDist-(tempDist*distTolerance))) 
            if (ScaledDist <= (tempDist+(tempDist*distTolerance)))

                %Compare angles only if not looking for first match
                if numTrue ~= 0

                    %Get angle between pair of stars in the input image
                    xdiff2 = centerCompare(1,1) - centerCompare(2,1);
                    ydiff2 = centerCompare(1,2) - centerCompare(2,2);
                    inputAngled = atan2d(ydiff2,xdiff2);

                    %Compare angles between template and input image
                    angleDiff = abs(inputAngled - templateAngled);

                    if angleDiff <= angleTolerance

                        %Set retArray regardless if this is the first match
                        if retArray(1) == 0
                            retArray = centerCompare; 
                            bestDistMatch = abs(ScaledDist - tempDist);
                            bestScaleMatch = abs(currSCF - tempSCF);
                            bestAngleMatch = angleDiff;

                        else

                            %Compare to previous best match
                            currDistMatch = abs(ScaledDist - tempDist);
                            currScaleMatch = abs(currSCF - tempSCF);
                            currAngleMatch = angleDiff;

                            distImprov = bestDistMatch/currDistMatch;
                            scaleImprov = bestScaleMatch/currScaleMatch;
                            angleImprov = bestAngleMatch/currAngleMatch;

                            score = distImprov + scaleImprov + angleImprov;
                            
                            % If this star is a better overall match
                            if score > 3 

                                %Assign values to retArray to be returned
                                retArray = centerCompare; 

                                %Set this star to the current best match
                                bestDistMatch = abs(ScaledDist - tempDist);
                                bestScaleMatch = abs(currSCF - tempSCF);
                                bestAngleMatch = angleDiff;

                            end
                        end

                    end

                %Skip angle comparison on the first match
                else

                    if retArray(1) == 0

                        %Assign values to retArray to be returned
                        retArray = centerCompare;

                        %Save info about current best match
                        bestDistMatch = abs(ScaledDist - tempDist);
                        bestScaleMatch = abs(currSCF - tempSCF);

                    else

                        %Compare to previous best match
                        currDistMatch = abs(ScaledDist - tempDist);
                        currScaleMatch = abs(currSCF - tempSCF);

                        distImprov = bestDistMatch/currDistMatch;
                        scaleImprov = bestScaleMatch/currScaleMatch;

                        score = distImprov + scaleImprov;

                        if score > 2 

                            %Assign values to retArray to be returned
                            retArray = centerCompare; 

                            %Set this star to the current best match
                            bestDistMatch = abs(ScaledDist - tempDist);
                            bestScaleMatch = abs(currSCF - tempSCF);

                        end

                    end

                end

            end
            end

        end
        end
    end
end