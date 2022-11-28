%%
pixel2mm('I:\Free-moving Feeding Data\20210323_AngelHair_PlexinD1CreER_Gcamp10\exp001', 'F:\Camera Parameters');
pixel2mm('I:\Free-moving Feeding Data\20210321_AngelHair_PlexinD1CreER_Gcamp11\exp001', 'F:\Camera Parameters');
pixel2mm('I:\Free-moving Feeding Data\20210323_AngelHair_PlexinD1CreER_Gcamp11\exp001', 'F:\Camera Parameters');
pixel2mm('I:\Free-moving Feeding Data\20210321_AngelHair_PlexinD1CreER_Gcamp12\exp001', 'F:\Camera Parameters');
pixel2mm('I:\Free-moving Feeding Data\20210323_AngelHair_PlexinD1CreER_Gcamp12\exp001', 'F:\Camera Parameters');
pixel2mm('D:\Free-moving Feeding Data\20201006_AngelHair_PlexinD1CreER_Gcamp7\exp001', 'F:\Camera Parameters');
pixel2mm('D:\Free-moving Feeding Data\20201006_AngelHair_PlexinD1CreER_Gcamp8\exp001', 'F:\Camera Parameters');
pixel2mm('D:\Free-moving Feeding Data\20200924_AngelHair_Fezf2CreER_Gcamp7\exp001', 'F:\Camera Parameters');

%%
pixel2mm('I:\Free-moving Feeding Data\20221009_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R1\exp001', 'F:\Camera Parameters @duke');
pixel2mm('I:\Free-moving Feeding Data\20221009_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R2\exp001', 'F:\Camera Parameters @duke');

%%
Raw_FPData_Preprocessing('D:\Fiber Photometry Data\20200801_AngelHair_Fezf2flpPlexinD1CreER_Gcamp4\exp001\signal0.csv',...
                         'D:\Fiber Photometry Data\20200801_AngelHair_Fezf2flpPlexinD1CreER_Gcamp4\exp001\trigger0.csv',...
                         'D:\Free-moving Feeding Data\20200801_AngelHair_Fezf2flpPlexinD1CreER_Gcamp4', [], 1, 1e13, 0, 10.6, 0, 1e10); % default are [], 1, 1e12, 0, 10.6, 0, 1e8

%%
Raw_FPData_Preprocessing('D:\Fiber Photometry Data\20221009_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R2\exp001\signal0.csv',...
                         'D:\Fiber Photometry Data\20221009_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R2\exp001\trigger0.csv',...
                         'I:\Free-moving Feeding Data\20221009_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R2', [0 1485; 2156 inf], 1, 1e13, 0, 10.6, 1, 1e10); % default are [], 1, 1e12, 0, 10.6, 0, 1e8

%%
Raw_FPData_Preprocessing('D:\Fiber Photometry Data\20221009_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R2\exp001\signal0.csv',...
                         'D:\Fiber Photometry Data\20221009_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R2\exp001\trigger0.csv',...
                         'I:\Free-moving Feeding Data\20221009_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R2', [], 1, 1e13, 0, 10.6, 1, 1e10); % default are [], 1, 1e12, 0, 10.6, 0, 1e8                    