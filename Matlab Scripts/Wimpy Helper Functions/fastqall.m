function [q_score, lengths, seq] = fastqall(path, prefix)
%Read in all fastq files in a directory and concatenate them into a single
%cell array

%INPUTS
%path - The directory where the fastq files are contained

%prefix - A filter used to identify fastq files from any other files
%present in the directory. Most nanopore outputs have "FAP" or "FAR" as a
%prefix

%OUTPUTS
%q_score - Cell array that stores the q score for each base for all reads
%length - vector containing lengths of all the reads
%seq - Cell array containing all the reads in the directory from all fastq
%files

%Get files in directory (specified by the "path" variable)
a = dir(path);

%If directory is not specified, set it to the current working directory, and
%throw a warning specifying that
if isempty(a)
    a = dir;
    warning(strcat(['Path not found: Using ', a(1).folder, ' as the path']))
end

%Initialize the size vector
sizes = zeros(length(a), 1);

for i = 1:length(a)
    sizes(i) = a(i).bytes;
end

seq = {};
q_score = {};

for i = 1:length(a) %for all the files in the directory
    n = a(i).name; %Store the name of the file
    if a(i).bytes > 10000 && sum(n(1:length(prefix)) == prefix) == length(prefix) %If the file size is > 10000 bytes and the prefix matches the prefix specified
        [~, s, q] = fastqread(n); %Read in the fastq
        seq(end+1:end+length(s), 1) = s; %Append to the seq variable
        q_score(end+1:end+length(q), 1) = q; %Append Q scores
    end
end

%Store lengths of all the reads
lengths = zeros(size(seq));
for i = 1:length(seq)
    lengths(i) = length(cell2mat(seq(i)));
end

end
