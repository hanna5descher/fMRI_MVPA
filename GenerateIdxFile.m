%inputs:
%fileName: GM mask file name (e.g., GM segmentation of a template). Must be
%.img file, each voxel is encoded as a unsigned 8-bit integer
%dim: x, y, and z dimensions of GM mask file in voxels, dim = [53 63 46]
%saveName: name of your idx file
%threshold, voxel values above which will be treated as GM
%r: radius of searchlight, in voxels
function saveData = GenerateIdxFile(fileName, dim, saveName, threshold, r)

fid = fopen(fileName, 'r');
data = fread(fid, inf, 'uint8');
fclose(fid);

data = data';
count = 0;
idx = zeros(size(data));
for i = 1 : length(data)
    if data(i) >= threshold
        count = count + 1;
        idx(i) = count;
    end
end

saveData = [count find(idx > 0)];
idx = reshape(idx, dim);

for i = 1 : dim(1)
    for j = 1 : dim(2)
        for k = 1 : dim(3)
            if idx(i, j, k) > 0
                %saveData(end + 1) = idx(i, j, k);
                %tIdx = [idx(i, j, k)];
                tIdx = [];
                for i0 = -r : r
                    for j0 = -r : r
                        for k0 = -r : r
                            if (i + i0 >= 1 && i + i0 <= dim(1) && j + j0 >= 1 && j + j0 <= dim(2) && k + k0 >= 1 && k + k0 <= dim(3) && i0 * i0 + j0 * j0 + k0 * k0 <= r * r)
                                if (idx(i + i0, j + j0, k + k0) > 0)
                                    tIdx(end + 1) = idx(i + i0, j + j0, k + k0);
                                end
                            end
                        end
                    end
                end
                tIdx = [length(tIdx) idx(i, j, k) tIdx];
                saveData = [saveData tIdx];
            end
        end
    end
end
        
fid = fopen(saveName, 'w');
fwrite(fid, saveData, 'int32');
fclose(fid);


