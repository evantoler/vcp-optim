function write_to_file(vec, fname)
% vec = vector of numerical data
% fname = name of file to write to

fid = fopen(fname, 'w'); % 'w' for write
for i = 1:numel(vec)
    fprintf(fid, '%1.15e\n', vec(i)); % print to double precision
end
fclose(fid);