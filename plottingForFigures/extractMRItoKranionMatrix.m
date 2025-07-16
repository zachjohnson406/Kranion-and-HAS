function [toKranion] = extractMRItoKranionMatrix(KRXfile,KRXpath)
KRXsubfolder = KRXfile(1:end-4);


disp('Parsing KRX data.')
%unzip([KRXpath KRXfile], KRXsubfolder);
cd([KRXpath KRXsubfolder]);

docNode = xmlread('model.xml');  % Reads the XML file

mrImages = docNode.getElementsByTagName('MR__Image');
toKranion = struct;
% Iterate through each <MR__Image> (assuming you may have multiple)
for i = 0:mrImages.getLength-1
    mrImage = mrImages.item(i);
    
    % Step 3: Find the <Attribute> node with key="ImageTransformMatExt"
    attributes = mrImage.getElementsByTagName('Attribute');
    for j = 0:attributes.getLength-1
        attribute = attributes.item(j);
        
        % Check if the 'key' attribute is 'ImageTransformMatExt'
        key = char(attribute.getAttribute('key'));
        if strcmp(key, 'ImageTransformMatExt')
            % Step 4: Extract all matrix elements dynamically
            matrixValues = NaN(4, 4);  % Initialize an empty 4x4 matrix
            
            % Loop through the elements m00 to m33
            for row = 0:3
                for col = 0:3
                     elementName = sprintf('m%01d%01d', row, col); % Creates names like 'm00', 'm01', ..., 'm33'
                    elementNode = attribute.getElementsByTagName(elementName).item(0);  % Get the node for that element
                    
                    % Extract and convert the value to double
                    if ~isempty(elementNode)
                        matrixValues(row+1, col+1) = str2double(char(elementNode.getTextContent()));
                    end
                end
            end
            
            % Display or store the matrix
            %disp('Extracted Matrix:');
            %disp(matrixValues);
            toKranion.(['matrix' num2str(i+1)]) = matrixValues;
        end
    end
end
end