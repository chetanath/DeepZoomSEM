function [status,Ztotal] = writeEDAXMultifieldMaps(outputfile,magnification,WD,x_vals,y_vals,z,stagetiltdegrees,stagerotationdegrees,amptime,dwelltime_microseconds,ev_per_channel,mapmatrix,MapPresetDwell,nLines,nPoints,nskip,CurrentWD,CurrentZ)

Z_v=str2double(z);
CurrentWD_v=str2double(CurrentWD);

Ztotal=[];

% if nargin<16
%     nskip=0;
% end

npts=numel(x_vals);

fid=fopen(outputfile,'w');
fprintf(fid,'<EDAX_TEAM_Multifield_Data_Table>\n');
fprintf(fid,'  <xs:schema id="EDAX_TEAM_Multifield_Data_Table" xmlns="" xmlns:xs="http://www.w3.org/2001/XMLSchema" xmlns:msdata="urn:schemas-microsoft-com:xml-msdata">\n');
fprintf(fid,'    <xs:element name="EDAX_TEAM_Multifield_Data_Table" msdata:IsDataSet="true" msdata:UseCurrentLocale="true"> \n');
fprintf(fid,'      <xs:complexType>\n');
fprintf(fid,'        <xs:choice minOccurs="0" maxOccurs="unbounded">\n');
fprintf(fid,'          <xs:element name="Version">\n');
fprintf(fid,'            <xs:complexType>\n');
fprintf(fid,'              <xs:sequence>\n');
fprintf(fid,'                <xs:element name="Version" type="xs:double" minOccurs="0" />\n');
fprintf(fid,'              </xs:sequence>\n');
fprintf(fid,'            </xs:complexType>\n');
fprintf(fid,'          </xs:element>\n');
fprintf(fid,'          <xs:element name="MultiFieldSettings">\n');
fprintf(fid,'            <xs:complexType>\n');
fprintf(fid,'              <xs:sequence>\n');
fprintf(fid,'                <xs:element name="Id" type="xs:string" />\n');
fprintf(fid,'                <xs:element name="StageOrdinal" type="xs:string" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="Type" type="xs:int" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="SubType" type="xs:int" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="SamplingRegionId" msdata:DataType="System.Guid, mscorlib, Version=2.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089" type="xs:string" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="Tag" msdata:DataType="System.Guid, mscorlib, Version=2.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089" type="xs:string" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="Text" type="xs:string" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="Magnification" type="xs:string" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="WorkingDistance" type="xs:float" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="StageX" type="xs:float" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="StageY" type="xs:float" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="StageZ" type="xs:float" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="StageTilt" type="xs:float" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="StageRotation" type="xs:float" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="AmpTime" type="xs:double" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="DwellTime" type="xs:double" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="EvPerChan" type="xs:double" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="MapMatrix" type="xs:int" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="PresetMode" type="xs:int" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="LineScan_ViewSizeX" type="xs:double" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="LineScan_ViewSizeY" type="xs:double" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="LineScan_X1" type="xs:double" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="LineScan_X2" type="xs:double" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="LineScan_Y1" type="xs:double" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="LineScan_Y2" type="xs:double" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="LineScan_Step" type="xs:double" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="LineScan_NPoints" type="xs:int" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="LineScan_FrameCount" type="xs:int" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="MapPresetDwell" type="xs:double" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="MapSize_Width" type="xs:double" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="MapSize_Height" type="xs:double" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="MapStartPoint_X" type="xs:double" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="MapStartPoint_Y" type="xs:double" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="Map_NLines" type="xs:int" minOccurs="0" />\n');
fprintf(fid,'                <xs:element name="Map_NPoints" type="xs:int" minOccurs="0" />\n');
fprintf(fid,'              </xs:sequence>\n');
fprintf(fid,'            </xs:complexType>\n');
fprintf(fid,'          </xs:element>\n');
fprintf(fid,'        </xs:choice>\n');
fprintf(fid,'      </xs:complexType>\n');
fprintf(fid,'      <xs:unique name="Constraint1" msdata:PrimaryKey="true">\n');
fprintf(fid,'        <xs:selector xpath=".//MultiFieldSettings" />\n');
fprintf(fid,'        <xs:field xpath="Id" />\n');
fprintf(fid,'      </xs:unique>\n');
fprintf(fid,'    </xs:element>\n');
fprintf(fid,'  </xs:schema>\n');
fprintf(fid,'  <Version>\n');
fprintf(fid,'    <Version>1</Version>\n');
fprintf(fid,'  </Version>\n');


for ii = 1:npts
    
    
    % Skip the first nskip maps
    if ii>nskip
        
        newZ=Z_v - (1e3*WD(ii)-CurrentWD_v); % Move Z stage instead
        Ztotal(end+1)=newZ;
        fprintf(fid,'  <MultiFieldSettings>\n');
        fprintf(fid,'    <Id>%d</Id>\n',ii);
        fprintf(fid,'    <StageOrdinal>%d</StageOrdinal>\n',ii);
        fprintf(fid,'    <Type>1</Type>\n');
        fprintf(fid,'    <SubType>4</SubType>\n');
        fprintf(fid,'    <SamplingRegionId>b3615ff6-d273-420b-99e7-ba758ce015ef</SamplingRegionId>\n');
        fprintf(fid,'    <Tag>46ec4958-f2a6-4813-9627-62865a08efbe</Tag>\n');
        fprintf(fid,'    <Magnification>%d</Magnification>\n',round(magnification));
        fprintf(fid,'    <WorkingDistance>%s</WorkingDistance>\n',CurrentWD); %1e3*WD(ii));
        fprintf(fid,'    <StageX>%.3f</StageX>\n',1e3*x_vals(ii));
        fprintf(fid,'    <StageY>%.3f</StageY>\n',1e3*y_vals(ii));
        fprintf(fid,'    <StageZ>%.3f</StageZ>\n',newZ);
        fprintf(fid,'    <StageTilt>%s</StageTilt>\n',stagetiltdegrees);
        fprintf(fid,'    <StageRotation>%s</StageRotation>\n',stagerotationdegrees);
        fprintf(fid,'    <AmpTime>%.2f</AmpTime>\n',amptime);
        fprintf(fid,'    <DwellTime>%.0E</DwellTime>\n',dwelltime_microseconds*1e-6);
        %        fprintf(fid,'    <Frames>%d</Frames>\n',framecount); % Doesn't work?? Must be set in TEAM...
        fprintf(fid,'    <EvPerChan>%d</EvPerChan>\n',ev_per_channel);
        fprintf(fid,'    <MapMatrix>6%d</MapMatrix>\n',mapmatrix);
        fprintf(fid,'    <PresetMode>1</PresetMode>\n');
        fprintf(fid,'    <MapPresetDwell>%.0E</MapPresetDwell>\n',MapPresetDwell);
        fprintf(fid,'    <MapSize_Width>1</MapSize_Width>\n');
        fprintf(fid,'    <MapSize_Height>1</MapSize_Height>\n');
        fprintf(fid,'    <MapStartPoint_X>0</MapStartPoint_X>\n');
        fprintf(fid,'    <MapStartPoint_Y>0</MapStartPoint_Y>\n');
        fprintf(fid,'    <Map_NLines>%d</Map_NLines>\n',round(nLines));
        fprintf(fid,'    <Map_NPoints>%d</Map_NPoints>\n',round(nPoints));
        fprintf(fid,'  </MultiFieldSettings>\n');
        
    end
end

fprintf(fid,'</EDAX_TEAM_Multifield_Data_Table>');

fprintf('Number of Regions: %d\n',npts)
fprintf('File written to: %s\n',outputfile)

status = fclose(fid);
