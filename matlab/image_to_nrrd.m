function parameters = image_to_nrrd(ct_img_path, nodule_img_path, dicom_path, pid)
    %% input part
    [lung_img_3d, nodule_img_3d, dicom_tags, thick, pixelsize, nodule_info] = fn_dicom_read(dicom_path,pid);
    fprintf('dicom images loaded ... \t\t\t %6.2f sec\n', toc);

    if numel(nodule_info) == 0
        return
    end
    [nodule_img_3d, nodule_info] = fn_nodule_info_update(lung_img_3d,nodule_img_3d,nodule_info,thick,pixelsize);

    %%
    meta = struct();
    meta.type = 'int16';
    meta.encoding = 'gzip';
    meta.spaceorigin = dicom_tags{1}.ImagePositionPatient';
    meta.spacedirections = [reshape(dicom_tags{1}.ImageOrientationPatient,3,2),[0;0;1]]*diag([pixelsize; thick]);
    meta.endian = 'little';

    parameters = table();
    parameters.PID = {pid};
    try
        parameters.PatientAge = dicom_tags{1}.PatientAge;
    catch
        parameters.PatientAge = '000Y';
    end
    try
        parameters.PatientSex = dicom_tags{1}.PatientSex;
    catch
        parameters.PatientSex = ' ';
    end
    parameters.KVP = dicom_tags{1}.KVP;
    parameters.XrayTubeCurrent = dicom_tags{1}.XRayTubeCurrent;
    parameters.Exposure = dicom_tags{1}.Exposure;
    try
        parameters.ExposureTime = dicom_tags{1}.ExposureTime;
    catch
        parameters.ExposureTime = 0;
    end
    try
        parameters.GeneratorPower = dicom_tags{1}.GeneratorPower;
    catch
        parameters.GeneratorPower = 0;
    end
    parameters.pixelsize = pixelsize';
    parameters.thick = thick;
    parameters.SliceThickness = dicom_tags{1}.SliceThickness;
    try
        parameters.SingleCollimationWidth = dicom_tags{1}.SingleCollimationWidth;
    catch
        parameters.SingleCollimationWidth = 0;
    end
    parameters.ConvolutionKernel = {dicom_tags{1}.ConvolutionKernel};
    try
        parameters.FilterType = {dicom_tags{1}.FilterType};
    catch
        parameters.FilterType = {''};
    end
    try
        parameters.ContrastBolusAgent = {dicom_tags{1}.ContrastBolusAgent};
        parameters.ContrastBolusRoute = {dicom_tags{1}.ContrastBolusRoute};
    catch
        parameters.ContrastBolusAgent = {''};
        parameters.ContrastBolusRoute = {''};
    end
    try
        parameters.StudyDescription = {dicom_tags{1}.StudyDescription};
    catch
        parameters.StudyDescription = {''};
    end
    try
        parameters.SeriesDescription = {dicom_tags{1}.SeriesDescription};
    catch
        parameters.SeriesDescription = {''};
    end
    parameters.Manufacturer = {dicom_tags{1}.Manufacturer};
    parameters.ManufacturerModelName = {dicom_tags{1}.ManufacturerModelName};
    try
        parameters.ProtocolName = {dicom_tags{1}.ProtocolName};
    catch
        parameters.ProtocolName = {''};
    end


    tic
    if ~isdir([nodule_img_path '/' pid]); mkdir([nodule_img_path '/' pid]); end
    fn_nrrdwrite([ct_img_path '/' pid  '_CT.nrrd'], int16(lung_img_3d(:,:,end:-1:1)), meta);
    writetable([nodule_info(:,[1:3 5:6 end 10:12 15:17]) nodule_info.Characteristics], [nodule_img_path '/' pid '/' pid '.csv']);

    meta.type = 'uint8';
    for sid = 1:4
        str_sid = num2str(sid);
        sid_nodule_image_3d =  uint8(bitand(uint8(nodule_img_3d),2^(sid-1))>0);
        sid_nodules = strcmp(nodule_info.sid, str_sid);
        if sum(sid_nodules) > 0
            fn_nrrdwrite([nodule_img_path '/' pid '/' pid '_CT_Phy' str_sid '-label.nrrd'], sid_nodule_image_3d(:,:,end:-1:1), meta);
        end
    end


    fprintf('dicom images and annotations converted ... \t\t\t %6.2f sec\n', toc);
end