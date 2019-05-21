function [ clust_indx ] = MS_find_clusters( tmp_ind, difvec )
%MS_FIND_CLUSTERS You will make me happy never reading this code!

clus_num = 1;

if isempty(difvec),
    clust_indx = [];  
else
 
    for dr=1:numel(difvec),
        if numel(difvec) == 1,
            if difvec == 1,
                clust_indx(clus_num) = {[tmp_ind(dr) tmp_ind(dr+1)]};  
            else
                clust_indx(clus_num) = {tmp_ind(dr)};
                clust_indx(clus_num+1) = {tmp_ind(dr+1)};
            end
            continue
        end        
        
        if (dr > 1 & difvec(dr-1) ~= 1 & difvec(dr) == 1) | (dr == 1 & difvec(dr) == 1),
            begin_val = tmp_ind(dr);
            if dr ~= numel(difvec),
                continue 
            end
        end

        if dr == numel(difvec),
            if difvec(dr) ~= 1,
                if difvec(dr-1) == 1,
                    end_val = tmp_ind(dr);
                    clust_indx(clus_num) = {[begin_val end_val]};
                    clust_indx(clus_num+1) = {tmp_ind(dr+1)};
                elseif difvec(dr-1) ~= 1,
                    clust_indx(clus_num) = {tmp_ind(dr+1)};
                end
            elseif difvec(dr) == 1,
                end_val = tmp_ind(dr+1);
                clust_indx(clus_num) = {[begin_val end_val]};
            end
            continue
        end

        if dr > 1 & difvec(dr-1) == 1 & difvec(dr) ~= 1,
            end_val = tmp_ind(dr);
            clust_indx(clus_num) = {[begin_val end_val]};
            clus_num = clus_num + 1;
            continue
        end

       if difvec(dr) ~= 1,
            clust_indx(clus_num) = {tmp_ind(dr)};
            clus_num = clus_num + 1;
            continue
       end    
    end

end

end

