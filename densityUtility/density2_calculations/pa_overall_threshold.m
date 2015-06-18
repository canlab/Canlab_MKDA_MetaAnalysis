
function pa_overall_threshold(studynames,allconditions,xyz,radius_mm,normby,mask)
%
% gets uncorrected and corrected p-values for prob of activation
%
% [allconditions,pointcounts,studynames,allcondindic,allcondnames] =
% dbcluster2indic(DB,DB,{testfield});
%
%

for i = 1:length(studynames)
    
        wh = find(strcmp(study,studynames{i}) & strcmp(allcond,allconditions{i}));
        XYZmm = xyz(wh,:);
        XYZmm(any(isnan(XYZmm),2),:) = [];
    
        if isempty(XYZmm)
            conmask = mask;
        else
            XYZvox = mm2vox(XYZmm,V.mat);   % get voxel coordinates
    
            conmask = xyz2mask(mask,XYZvox);   % put points in mask - could weight by Z-scores here, if desired.
            spm_smooth(conmask,conmask,radius_mm); % smooth it!
            
            conmask = conmask ./ normby;    % normalize so center of study activation = 1
            conmask(conmask > 1) = 1;       % limit to max of 1, in case multiple nearby points in same contrast
                                            % max activation for a single
                                            % contrast is 1.
    
        end
        
        str = [studynames{i} '_' testfield '_' allconditions{i} '.img'];
        V.fname = str;
        
        V.descrip = allconditions{i};
        %str = ['V.' testfield ' = ''' allconditions{i} ''';'];
        %eval(str)
        
        
        warning off, spm_write_vol(V,conmask);, warning on
        PP = str2mat(PP,str);
    end % loop thru contrasts
    PP = PP(2:end,:);
else
    % we have filenames already created
end

fprintf(1,'Done %3.0f images in %3.0f s',length(studynames),etime(clock,t1))
    
OUT.PP = PP;