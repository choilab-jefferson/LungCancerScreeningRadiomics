%% set environment
% external toolbox path
addpath(genpath('toolbox'))


%% lung image analysis
% lung image analysis path
LIA_path = 'toolbox/lung-image-analysis';

% set global values
global dicom_path;
global data_path;
global output_path;

dicom_path = '../DATA/LIDC-IDRI'; %dcm files directory
data_path = '../DATA/'; 
output_path = 'output/';

ConformalizedMCF = ['docker run -v ' tempdir ':' tempdir ' wookjinchoi/conformalized_mcf:latest ConformalizedMCF'];

% 
esph_factor = (3/(4*pi))^(1/3);
iso = ''; % '-iso'; '-1mm';


% module path values
util_path=[LIA_path '/util'];
input_path=[LIA_path '/io'];
interpolation_path=[LIA_path '/interpolation'];
segmentation_path=[LIA_path '/lung_segmentation'];
nodule_segmentation_path=[LIA_path '/nodule_segmentation'];
candidates_path=[LIA_path '/nodule_candidate_detection'];
feature_extraction_path=[LIA_path '/feature_extraction'];
evaluation_path=[LIA_path '/evaluation'];

% module addpath
addpath(genpath(util_path));
addpath(genpath(input_path));
addpath(genpath(interpolation_path));
addpath(genpath(segmentation_path));
addpath(genpath(nodule_segmentation_path));
addpath(genpath(candidates_path));
addpath(genpath(feature_extraction_path));
addpath(genpath(evaluation_path));


%% selected patients who have ground truth
% primary lung cancer, metastasis and bvoxelize_meshesenign
selected_157 = {'LIDC-IDRI-0068'
    'LIDC-IDRI-0071'
    'LIDC-IDRI-0072'
    'LIDC-IDRI-0088'
    'LIDC-IDRI-0090'
    'LIDC-IDRI-0091'
    'LIDC-IDRI-0100'
    'LIDC-IDRI-0118'
    'LIDC-IDRI-0124'
    'LIDC-IDRI-0129'
    'LIDC-IDRI-0135'
    'LIDC-IDRI-0137'
    'LIDC-IDRI-0138'
    'LIDC-IDRI-0143'
    'LIDC-IDRI-0149'
    'LIDC-IDRI-0159'
    'LIDC-IDRI-0161'
    'LIDC-IDRI-0162'
    'LIDC-IDRI-0163'
    'LIDC-IDRI-0164'
    'LIDC-IDRI-0165'
    'LIDC-IDRI-0166'
    'LIDC-IDRI-0167'
    'LIDC-IDRI-0168'
    'LIDC-IDRI-0169'
    'LIDC-IDRI-0171'
    'LIDC-IDRI-0173'
    'LIDC-IDRI-0174'
    'LIDC-IDRI-0175'
    'LIDC-IDRI-0176'
    'LIDC-IDRI-0178'
    'LIDC-IDRI-0179'
    'LIDC-IDRI-0180'
    'LIDC-IDRI-0181'
    'LIDC-IDRI-0182'
    'LIDC-IDRI-0183'
    'LIDC-IDRI-0184'
    'LIDC-IDRI-0185'
    'LIDC-IDRI-0186'
    'LIDC-IDRI-0187'
    'LIDC-IDRI-0188'
    'LIDC-IDRI-0189'
    'LIDC-IDRI-0190'
    'LIDC-IDRI-0191'
    'LIDC-IDRI-0192'
    'LIDC-IDRI-0193'
    'LIDC-IDRI-0194'
    'LIDC-IDRI-0197'
    'LIDC-IDRI-0198'
    'LIDC-IDRI-0200'
    'LIDC-IDRI-0202'
    'LIDC-IDRI-0203'
    'LIDC-IDRI-0205'
    'LIDC-IDRI-0207'
    'LIDC-IDRI-0210'
    'LIDC-IDRI-0211'
    'LIDC-IDRI-0212'
    'LIDC-IDRI-0213'
    'LIDC-IDRI-0214'
    'LIDC-IDRI-0217'
    'LIDC-IDRI-0220'
    'LIDC-IDRI-0221'
    'LIDC-IDRI-0222'
    'LIDC-IDRI-0223'
    'LIDC-IDRI-0224'
    'LIDC-IDRI-0225'
    'LIDC-IDRI-0226'
    'LIDC-IDRI-0230'
    'LIDC-IDRI-0231'
    'LIDC-IDRI-0232'
    'LIDC-IDRI-0233'
    'LIDC-IDRI-0234'
    'LIDC-IDRI-0235'
    'LIDC-IDRI-0236'
    'LIDC-IDRI-0237'
    'LIDC-IDRI-0239'
    'LIDC-IDRI-0242'
    'LIDC-IDRI-0243'
    'LIDC-IDRI-0244'
    'LIDC-IDRI-0245'
    'LIDC-IDRI-0246'
    'LIDC-IDRI-0247'
    'LIDC-IDRI-0248'
    'LIDC-IDRI-0249'
    'LIDC-IDRI-0250'
    'LIDC-IDRI-0251'
    'LIDC-IDRI-0252'
    'LIDC-IDRI-0253'
    'LIDC-IDRI-0254'
    'LIDC-IDRI-0255'
    'LIDC-IDRI-0256'
    'LIDC-IDRI-0257'
    'LIDC-IDRI-0258'
    'LIDC-IDRI-0260'
    'LIDC-IDRI-0261'
    'LIDC-IDRI-0264'
    'LIDC-IDRI-0265'
    'LIDC-IDRI-0266'
    'LIDC-IDRI-0267'
    'LIDC-IDRI-0268'
    'LIDC-IDRI-0270'
    'LIDC-IDRI-0271'
    'LIDC-IDRI-0272'
    'LIDC-IDRI-0273'
    'LIDC-IDRI-0274'
    'LIDC-IDRI-0275'
    'LIDC-IDRI-0276'
    'LIDC-IDRI-0277'
    'LIDC-IDRI-0278'
    'LIDC-IDRI-0279'
    'LIDC-IDRI-0280'
    'LIDC-IDRI-0281'
    'LIDC-IDRI-0282'
    'LIDC-IDRI-0283'
    'LIDC-IDRI-0285'
    'LIDC-IDRI-0286'
    'LIDC-IDRI-0287'
    'LIDC-IDRI-0288'
    'LIDC-IDRI-0289'
    'LIDC-IDRI-0290'
    'LIDC-IDRI-0314'
    'LIDC-IDRI-0325'
    'LIDC-IDRI-0332'
    'LIDC-IDRI-0377'
    'LIDC-IDRI-0385'
    'LIDC-IDRI-0399'
    'LIDC-IDRI-0405'
    'LIDC-IDRI-0454'
    'LIDC-IDRI-0470'
    'LIDC-IDRI-0493'
    'LIDC-IDRI-0510'
    'LIDC-IDRI-0522'
    'LIDC-IDRI-0543'
    'LIDC-IDRI-0559'
    'LIDC-IDRI-0562'
    'LIDC-IDRI-0568'
    'LIDC-IDRI-0576'
    'LIDC-IDRI-0580'
    'LIDC-IDRI-0610'
    'LIDC-IDRI-0624'
    'LIDC-IDRI-0766'
    'LIDC-IDRI-0771'
    'LIDC-IDRI-0772'
    'LIDC-IDRI-0811'
    'LIDC-IDRI-0818'
    'LIDC-IDRI-0875'
    'LIDC-IDRI-0893'
    'LIDC-IDRI-0905'
    'LIDC-IDRI-0921'
    'LIDC-IDRI-0924'
    'LIDC-IDRI-0939'
    'LIDC-IDRI-0965'
    'LIDC-IDRI-0994'
    'LIDC-IDRI-1002'
    'LIDC-IDRI-1004'
    'LIDC-IDRI-1010'
    'LIDC-IDRI-1011'};

% primary lung cancer and benign
selected_72 = {'LIDC-IDRI-0072'
    'LIDC-IDRI-0090'
    'LIDC-IDRI-0138'
    'LIDC-IDRI-0149'
    'LIDC-IDRI-0162'
    'LIDC-IDRI-0163'
    'LIDC-IDRI-0166'
    'LIDC-IDRI-0167'
    'LIDC-IDRI-0168'
    'LIDC-IDRI-0171'
    'LIDC-IDRI-0178'
    'LIDC-IDRI-0180'
    'LIDC-IDRI-0183'
    'LIDC-IDRI-0185'
    'LIDC-IDRI-0186'
    'LIDC-IDRI-0187'
    'LIDC-IDRI-0191'
    'LIDC-IDRI-0203'
    'LIDC-IDRI-0211'
    'LIDC-IDRI-0212'
    'LIDC-IDRI-0233'
    'LIDC-IDRI-0234'
    'LIDC-IDRI-0242'
    'LIDC-IDRI-0246'
    'LIDC-IDRI-0247'
    'LIDC-IDRI-0249'
    'LIDC-IDRI-0256'
    'LIDC-IDRI-0257'
    'LIDC-IDRI-0265'
    'LIDC-IDRI-0267'
    'LIDC-IDRI-0268'
    'LIDC-IDRI-0270'
    'LIDC-IDRI-0271'
    'LIDC-IDRI-0273'
    'LIDC-IDRI-0275'
    'LIDC-IDRI-0276'
    'LIDC-IDRI-0277'
    'LIDC-IDRI-0283'
    'LIDC-IDRI-0286'
    'LIDC-IDRI-0289'
    'LIDC-IDRI-0290'
    'LIDC-IDRI-0314'
    'LIDC-IDRI-0325'
    'LIDC-IDRI-0332'
    'LIDC-IDRI-0377'
    'LIDC-IDRI-0385'
    'LIDC-IDRI-0399'
    'LIDC-IDRI-0405'
    'LIDC-IDRI-0454'
    'LIDC-IDRI-0470'
    'LIDC-IDRI-0493'
    'LIDC-IDRI-0510'
    'LIDC-IDRI-0522'
    'LIDC-IDRI-0543'
    'LIDC-IDRI-0559'
    'LIDC-IDRI-0562'
    'LIDC-IDRI-0568'
    'LIDC-IDRI-0580'
    'LIDC-IDRI-0610'
    'LIDC-IDRI-0624'
    'LIDC-IDRI-0766'
    'LIDC-IDRI-0771'
    'LIDC-IDRI-0811'
    'LIDC-IDRI-0875'
    'LIDC-IDRI-0905'
    'LIDC-IDRI-0921'
    'LIDC-IDRI-0924'
    'LIDC-IDRI-0939'
    'LIDC-IDRI-0965'
    'LIDC-IDRI-0994'
    'LIDC-IDRI-1002'
    'LIDC-IDRI-1004'};


%% structural elements for morphology
[x,y,z] = ndgrid(-1:1);
se1 = strel(sqrt(x.^2 + y.^2 + z.^2) <=1);

[x,y,z] = ndgrid(-2:2);
se2 = strel(sqrt(x.^2 + y.^2 + z.^2) <=2);

[x,y,z] = ndgrid(-4:4);
se4 = strel(sqrt(x.^2 + y.^2 + z.^2) <=4);

