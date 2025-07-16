function sonicationData = extractSonicationData(xmlFile)
    % Parse the XML file into a structure
    xmlStruct = xml2struct(xmlFile);
    
    % Initialize an empty cell array to store the extracted data
    sonicationData = {};
    
    % Get all Sonication tags (assuming there are multiple Sonication tags)
    sonications = xmlStruct.Children;  % Top level of the XML
    
    % Loop over each Sonication tag
    for i = 1:length(sonications)
        if strcmp(sonications(i).Name, 'Sonication')
            % Extract the attributes from the Sonication tag
            num = str2double(sonications(i).Attributes.num);
            power_w = str2double(sonications(i).Attributes.power__w);
            duration_s = str2double(sonications(i).Attributes.duration__s);
            frequency = str2double(sonications(i).Attributes.frequency);
            
            % Find the ThermometryPhase tag within Sonication
            thermometryPhase = getElement(sonications(i), 'ThermometryPhase');
            
            % Find the Attributes tag within ThermometryPhase
            attributes = getElement(thermometryPhase, 'Attributes');
            
            % Find the Attribute tag with key 'ImagePosition'
            imagePositionAttr = getAttributeByKey(attributes, 'ImagePosition');
            
            % Extract the float values (there should be 3 <float> tags)
            floats = extractFloatValues(imagePositionAttr);
            
            % Store the data for this Sonication in the result cell array
            sonicationData{end+1} = struct('num', num, ...
                                           'power_w', power_w, ...
                                           'duration_s', duration_s, ...
                                           'frequency', frequency, ...
                                           'imagePosition', floats);
        end
    end
end

function element = getElement(parent, tagName)
    % Helper function to get the child element by tag name
    element = [];
    for j = 1:length(parent.Children)
        if strcmp(parent.Children(j).Name, tagName)
            element = parent.Children(j);
            break;
        end
    end
end

function attribute = getAttributeByKey(attributes, key)
    % Helper function to find an attribute by its 'key' attribute
    attribute = [];
    for k = 1:length(attributes.Children)
        if isfield(attributes.Children(k).Attributes, 'key') && ...
           strcmp(attributes.Children(k).Attributes.key, key)
            attribute = attributes.Children(k);
            break;
        end
    end
end

function floatValues = extractFloatValues(attribute)
    % Extract the <float> values from an Attribute element
    floatValues = [];
    for m = 1:length(attribute.Children)
        if strcmp(attribute.Children(m).Name, 'float')
            floatValues = [floatValues; str2double(attribute.Children(m).Text)];
        end
    end
end
