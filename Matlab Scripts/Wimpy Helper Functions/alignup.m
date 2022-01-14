function [perc_R,perc_F, new_seq, right_seq, flip, start_points] = alignup(q,ref,thresh)
% Align reads and flip reads as necessary, by aligning the reads to a small
% reference sequence (~100bp), and flipping a read if it aligns well to the
% reverse complement of the ref seq. The threshold for considering an
% alignment as a "hit" is user defined by the variable "thresh"

%INPUTS
% q - Cell array with reads (1 x n or n x 1 array, where n is the number of reads)
% ref - reference sequence to align reads to
% thresh - alignment threshold for considering an alignment a hit
% (fraction, between 0 and 1, where 1 represents perfect alignment)

%OUTPUTS

%right_seq - Cell array with all reads that align to the reverse complement
% flipped (and others kept as is), so the top strand is being read from 5' to 3'

%flip - Vector that indicates whether a read was flipped or not - 
%   1 = aligned to reverse (flipped and stored in new_seq)
%   2 = aligned to forward (stored in new_seq as is)
%   -1 = does not align to reference (stored in new_seq as 'X')

%new_seq - Cell array similar to right_seq, with the exception that reads
%that do not meet the alignment threshold being replaced by 'X'. Further,
%the reads are re-formatted to start from the reference sequence (instead
%of their original start points)

%perc_R - Alignment fraction (between 0-1) to the reverse of the reference
%seq

%perc_F - Alignment fraction to reference, among reads that were not
%flipped (note: reads that aligned to the reverse complement will not be
%checked for alignment here, and will be stored as 0s in this variable

%Start_points - the index of the start point of the alignment to the reference, 
%for every read

%If the reference sequence is > 150 bp long, trim and keep only the first
%150 bp

if length(ref) > 150
    ref = upper(ref(1:150));
end

%Reverse complement of the reference sequence
ref2 = seqrcomplement(ref);

%Initialize all variables
flip = zeros(size(q)); %Initialize the flip vector
perc_R = zeros(size(q)); %sequences aligned to reverse complement of ref
perc_F = zeros(size(q)); %sequences aligned to forward strand of ref
right_seq = cell(size(q)); %sequence list with non-aligning sequences left in
start_points = zeros(size(q));
new_seq = cell(size(q)); %sequence list with non-aligning sequences replaced with 'X'

for i = 1:length(q) %For every read
    if mod(i, 1000) == 0 %Display progress in chunks of 1000 reads
        disp(i)
    end
    a = cell2mat(q(i)); %Read in sequence
    [~, x1, y1] = swalign(a, ref2, 'alphabet', 'NT'); %Align to REVERSE of reference
    perc_R(i) = sum(x1(2, :) == '|')/length(ref2); %compute alignment fraction
    
    if perc_R(i) > thresh %If alignment is greater than threshold
         flip(i) = 1; %Flip the read
         a = seqrcomplement(a); %Save reverse complement of read
         right_seq(i) = cellstr(a);
         start_points(i) = y1(1)+length(ref); %Save the start point
         
         d = length(a) - (y1(1)+length(ref)); %Distance from alignment point to end
         if d > 0
            a2 = a(d:end);
            a = strcat(a2, a(1:d-1)); %Re arrange the read so it starts from the alignment point
            new_seq(i) = cellstr(a);
         else
             new_seq(i) = cellstr('X');
         end
    end

     if flip(i) == 0 %If read was not flipped, it should align to the forward 
         [~, x2, y2] = swalign(a, ref, 'alphabet', 'NT'); %Align to forward
         perc_F(i) = sum(x2(2, :) == '|')/length(ref); %Check alignment percentage (only in unflipped reads)
         if perc_F(i) > thresh %If greater than threshold
             flip(i) = 2; %Mark that these are good reads
             right_seq(i) = cellstr(a); %Save as-is (no flipping)
             start_points(i) = y2(1); %Save start point
             
             ups = a(1:y2(1)); %From first base to reference start point
             downs = a(y2(1):end); %From ref start point to end
             a = strcat(downs, ups); %Re-arrange read
             new_seq(i) = cellstr(a);

        else
             right_seq(i) = cellstr(a);
             flip(i) = -1;
             start_points(i) = -1;
             new_seq(i) = cellstr('X');
         end
     
    end

end
end
