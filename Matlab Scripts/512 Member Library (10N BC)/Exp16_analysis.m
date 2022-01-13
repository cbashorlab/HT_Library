%% Analysing Experiment 16 data (pKR194 v3)

%Add the directory with WIMPY helper functions to the path
addpath('C:\Lab Users\ROC\Wimpy Helper Functions');

%load in data using fastqall
[~, l, seq] = fastqall('x', 'fastq');
histogram(l);

%Plot read length histogram
figure
histogram(l)
xlabel('Length (kb)');
xlim([0 20000]);
ylabel('Read Count');
%ylim([0 4500000]);

%filter reads by expected length
q = seq(l > 7000 & l < 14000); %for 512-member library

%flip and index reads using alignup function
ref = upper('agggccaccttcatgagctacaacaccatcatcagcaacagcctgagcttcgacatcgtgaacaagagcctgcagttcaagtacaagacccagaaggccaccatcctggaggccagcctgaagaagctgatccccgcctgggagttcaccat');
thresh = 0.75;
[~, ~, pKR194, ~, ~, ~] = alignup(q,ref,thresh);
pKR194 = pKR194(~contains(pKR194, 'X'));


%% Indexing reads to Puro

% ref = upper('cggccgaacgagcaggcgtcccagccttcttggaAacgagcgcacccaggaatcttccattctacgagcgactgggtttcactgtgacagcagacgtagaagtccctgaagggcctcggacttggtgcatgaccagaaagcctggggcttg');
% perc_x3 = zeros(size(pKR194));
% new_seq = cell(size(pKR194));
% start_points = zeros(size(pKR194));
% 
% for i = 1:length(pKR194)
%     
%     if mod(i, 500) == 0
%     disp(i)
%     end
%     
%     a = cell2mat(pKR194(i));
%     [~, x1, y1] = swalign(a, ref, 'alphabet', 'NT');
%     perc_x3(i) = sum(x1(2, :) == '|')/length(ref);
% 
%     start_points(i) = y1(1)+length(ref);
%     d = length(a) - (y1(1)+length(ref));
%     
%     if d > 0
%         a2 = a(d:end);
%         a = strcat(a2, a(1:d-1));
%         new_seq(i) = cellstr(a);
%     else
%         new_seq(i) = cellstr('X');
%     end
%     
% end

%% EU ASSIGNMENT

%Align reads to mCherry and index mCherry location
mCh_align_perc = zeros(size(pKR194));
mCh_align_start = mCh_align_perc;

mCh = upper('ggtgaggaagaTaacatggcgattataaaggaatttatgcggttcaaggtgcacatggagggttcagttaatggacacgagttcgaaatcgaaggtgagg');

for i = 1:length(pKR194)
    if mod(i, 1000) == 0
       disp(i)
    end
    a = cell2mat(pKR194(i));
    
    [~, x, y] = swalign(a, mCh, 'Alphabet', 'NT');
    mCh_align_perc(i) = sum(x(2, :) == '|')/length(mCh);
    mCh_align_start(i) = y(1);
end


mCh_align_start = mCh_align_start(mCh_align_start > 0);

%Check for Promoters and Kozaks in regions UPSTREAM of mCherry start site
%Align to all EU components using the tiling method
%Create upstream and downstream mCherry regions

p_region = cell(size(pKR194));
t_region = p_region;

for j = 1:length(pKR194)
    a = cell2mat(pKR194(j));
    if mCh_align_start(j) < 3000 && mCh_align_start(j) > 0
        p_region(j) = cellstr(a(1:mCh_align_start(j)));
        t_region(j) = cellstr(a(mCh_align_start(j):mCh_align_start(j) + 2000));
    else
         p_region(j) = cellstr('X');
         t_region(j) = cellstr('X');
    end
end

pKR194 = pKR194(~contains(p_region, 'X'));
p_region = p_region(~contains(p_region, 'X'));
t_region = t_region(~contains(t_region, 'X'));

%PROMOTERS (containment)
%Read in promoter reference files
promoters = readcell('512-Promoters.xlsx');
promoters = promoters(2:9, :);
promoters = promoters(:, 2);
promoter_reads = zeros(length(pKR194), 8); %Initialize matrix to store number of tiles contained for each promoter for every read

for p = 1:8
    disp(p)
    a = upper(cell2mat(promoters(p)));
    for m = 1:length(a) - 12  %Break into 12bp tiles
        b = a(m:m+11);
        promoter_reads(:, p) = promoter_reads(:, p) + contains(p_region, b);  %Add 1 for the reads that contained the tile
    end
end

%Scale the number of tiles found by the length of the promoters, so every
%read is normalized against it (and represents what fraction of the tiles
%were found for each promoter)

promoter_lengths = repmat([length(cell2mat(promoters(1))), length(cell2mat(promoters(2))), length(cell2mat(promoters(3))), length(cell2mat(promoters(4))), length(cell2mat(promoters(5))), length(cell2mat(promoters(6))), length(cell2mat(promoters(7))), length(cell2mat(promoters(8)))], length(pKR194), 1);
promoter_reads_scaled = promoter_reads./promoter_lengths;

%Confusion matrix for promoter assignment
pconf = zeros(8, 8);
for i = 1:8
    for j = 1:8
        pconf(i, j) = sum(promoter_reads_scaled(:, i) > 0.04 & promoter_reads_scaled(:, j) > 0.04);
    end
end

%Promoters 1 and 2 are very similar and thus have higher confusion.
%Increase the threshold for those two promoters

for i = 1:2
    for j = 1:2
        pconf(i, j) = sum(promoter_reads_scaled(:, i) > 0.3 & promoter_reads_scaled(:, j) > 0.3);
    end
end
figure
confusionchart(pconf)

%[p_score, max_p] = max(promoter_reads_scaled'); %re-run later


%TERMINATORS
terminators = readcell('512-Terminators.xlsx');
terminators = terminators(:, 2);
terminator_reads = zeros(length(t_region), 8);

for t = 1:8
    disp(t)
    a = upper(cell2mat(terminators(t)));
    for m = 1:length(a) - 12
        b = a(m:m+11);
        terminator_reads(:, t) = terminator_reads(:, t) + contains(t_region, b);
    end
end

terminator_lengths = repmat([length(cell2mat(terminators(1))), length(cell2mat(terminators(2))), length(cell2mat(terminators(3))), length(cell2mat(terminators(4))), length(cell2mat(terminators(5))), length(cell2mat(terminators(6))), length(cell2mat(terminators(7))), length(cell2mat(terminators(8)))], length(t_region), 1);
terminator_reads_scaled = terminator_reads./terminator_lengths;

tconf = zeros(8, 8);
for i = 1:8
    for j = 1:8
        tconf(i, j) = sum(terminator_reads_scaled(:, i) > 0.15 & terminator_reads_scaled(:, j) > 0.15);
    end
end

for i = 1:2
    for j = 1:2
        tconf(i, j) = sum(terminator_reads_scaled(:, i) > 0.30 & terminator_reads_scaled(:, j) > 0.30);
    end
end

figure
confusionchart(tconf)


%KOZAKS - look in p_region
test_kozaks = readcell('512-Kozak-Tests.xlsx');
test_kozaks = test_kozaks(:, 5);

%Generating confusion matrices
kozak_conf = zeros(8, 8);

for i = 1:8
    for j = 1:8
        kozak_conf(i, j) = sum(contains(p_region, cell2mat(test_kozaks(i))) & contains(p_region, cell2mat(test_kozaks(j))));
    end
end

figure
confusionchart(kozak_conf)

%Fix kozak 2 and kozak 5 confusion
test_kozaks2 = readcell('512-Kozak-Tests.xlsx');
test_kozaks2 = test_kozaks2(:, 6);

for i = [2, 5]
    for j = [2, 5]
        kozak_conf(i, j) = sum(contains(p_region, cell2mat(test_kozaks2(i))) & contains(p_region, cell2mat(test_kozaks2(j))));
    end
end    

figure
confusionchart(kozak_conf)

%Assigning Reads
kassign = zeros(length(pKR194), 1);

for i = 1:length(t_region)
    disp(i)
    a = cell2mat(p_region(i));
    for j = [1, 3, 4, 6:8]
        if contains(a, cell2mat(test_kozaks(j)))
%            pKR194_kassigned(i) = pKR194(i);
            kassign(i) = j;
            break
        end
    end
end

kassign(contains(p_region, test_kozaks(2)) & kassign == 0) = 2;
kassign(contains(p_region, test_kozaks(5)) & ~contains(p_region, test_kozaks(2)) & kassign == 0) = 5;


%TOTAL EU COMPOSITION
[p_score, max_p] = max(promoter_reads_scaled');
passign = max_p';

[t_score, max_t] = max(terminator_reads_scaled');
tassign = max_t';

%ASSIGN EUs in 8x8x8 matrix
EU_assignments = zeros(8, 8, 8);

for i = 1:8
    for j = 1:8
        for k = 1:8
            EU_assignments(i, j, k) = sum(passign == i & kassign == j & tassign == k);
        end
    end
end

%% Barcoding

barcodes = cell(size(pKR194));
both_scars = contains(pKR194, 'AACGTCGC') & contains(pKR194, 'GAATTC');
len_bc = zeros(size(pKR194)) - 1;

for i = 1:length(barcodes)
    disp(strcat('BC_', num2str(i)))
    if both_scars(i)
       a = cell2mat(pKR194(i));
       found = strfind(a, 'AACGTCGC');
       found = found(found > 2000); % & found < 4000);
       for j = 1:length(found)
           if length(a) < found+50
               x = a(found(j):end);
           else
           x = a(found(j):found(j)+50);
           found_3 = strfind(x, 'GAATT');
           if isempty(found_3) == 0
               if length(found_3) > 1
                   found_3 = found_3(1);
               end
               barcodes(i) = cellstr(x(9:found_3 - 1));
               len_bc(i) = length(x(9:found_3-1));
               break
           end
           end
        
        end
    end
end


bc_10bp = barcodes(tassign > 0 & passign > 0 & kassign > 0 & len_bc == 10); % & promoter_index > 0);
terminators_barcoded = tassign(len_bc == 10 & tassign > 0 & passign > 0 & kassign > 0);
promoters_barcoded = passign(len_bc == 10 & tassign > 0 & passign > 0 & kassign > 0);
kozaks_barcoded = kassign(len_bc == 10 & tassign > 0 & passign > 0 & kassign > 0);

%generate list of unique barcodes
[bc_unique, c, b] = unique(bc_10bp); %Find the sequences of unique RBSs (RBS_unique), and their indices (b)
bc_count = zeros(1, length(bc_unique));
eu_dist = zeros(512, length(bc_unique));
variant_barcoded = (promoters_barcoded - 1)*64 + (kozaks_barcoded - 1)*8 + terminators_barcoded;

for i = 1:length(bc_unique)
    bc_count(i) = sum(b == i);
    x = variant_barcoded(b == i);
    for j = 1:512
        eu_dist(j, i) = sum(x == j);
    end
end

figure
binsdetected = sum(eu_dist > 0);
histogram(binsdetected)

figure
for i = 1:5
    subplot(5, 1, i)
    histogram(binsdetected(bc_count > i-1))
end

%cross reference with barcodes from Illumina run
histogram(bc_count)

for j = 1:length(bc_unique)
    
   x(j) = sum(contains(Barcodes_filtered, cell2mat(bc_unique(j))));
   
end

bc_found = bc_unique(x > 0);
bc_found_count = bc_count(x > 0);

eu_dist_found = eu_dist(:, x > 0);

binsdetected_found = sum(eu_dist_found > 0);
histogram(binsdetected_found)

for i = 1:10
subplot(10, 1, i)
histogram(binsdetected_found(bc_found_count > i-1))

end


%% Barcode -> EU assignment without Kozaks

% bc_10bp_noK = barcodes(tassign > 0 & passign > 0 & len_bc == 10); %kassign > 0; % & promoter_index > 0);
% terminators_barcoded_noK = tassign(len_bc == 10 & tassign > 0 & passign > 0);
% promoters_barcoded_noK = passign(len_bc == 10 & tassign > 0 & passign > 0);
% 
% %generate list of unique barcodes
% [bc_unique_noK, c_noK, b_noK] = unique(bc_10bp_noK); %Find the sequences of unique RBSs (RBS_unique), and their indices (b)
% bc_count_noK = zeros(1, length(bc_unique_noK));
% eu_dist_noK = zeros(64, length(bc_unique_noK));
% variant_barcoded_noK = (promoters_barcoded_noK - 1)*8 + terminators_barcoded_noK;
% 
% for i = 1:length(bc_unique_noK)
%     bc_count_noK(i) = sum(b_noK == i);
%     x_noK = variant_barcoded_noK(b_noK == i);
%     for j = 1:64
%         eu_dist_noK(j, i) = sum(x_noK == j);
%     end
% end
% 
% 
% binsdetected_noK = sum(eu_dist_noK > 0);
% histogram(binsdetected_noK)
% 
% figure
% for i = 1:5
%     subplot(5, 1, i)
%     histogram(binsdetected_noK(bc_count_noK > i-1))
% end
% 
% %cross reference with barcodes from Illumina run
% histogram(bc_count_noK)
% 
% for j = 1:length(bc_unique_noK)
%     
%    x_noK(j) = sum(contains(Barcodes_filtered, cell2mat(bc_unique_noK(j))));
%    
% end
% 
% bc_found_noK = bc_unique_noK(x_noK > 0);
% bc_found_count_noK = bc_count_noK(x_noK > 0);
% 
% eu_dist_found_noK = eu_dist_noK(:, x_noK > 0);
% 
% binsdetected_found_noK = sum(eu_dist_found_noK > 0);
% figure
% histogram(binsdetected_found_noK)
% 
% for i = 1:10
% subplot(10, 1, i)
% histogram(binsdetected_found_noK(bc_found_count_noK > i-1))
% 
% end

%% PCA for Kozak Assignment

%Get 30bp upstream from mCh Start

% truncation = 30;
% mCh_start = mCh_align_start(mCh_align_start < 3000 & mCh_align_start > 0);
% reads_kozak_region = cell(size(mCh_start));
% 
% for i = 1:length(reads_kozak_region)
% %    if mod(i, 10000) == 0
% %        disp(i)
% %    end
%     a = cell2mat(pKR194(i));
%     if mCh_start(i) > 1500
%     reads_kozak_region(i) = cellstr(a(mCh_start(i) - 29:mCh_start(i)));
%     else
%         reads_kozak_region(i) = cellstr('X');
%     end
% end
%     
% reads_kozak_region = reads_kozak_region(~contains(reads_kozak_region, 'X'));







