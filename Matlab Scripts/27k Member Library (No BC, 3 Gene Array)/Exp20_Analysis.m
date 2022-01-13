%% Analysing Experiment 20 data (27k library)

%Add the directory with WIMPY helper functions to the path
addpath('C:\Lab Users\ROC\Wimpy Helper Functions');

%Load directory
direc = dir;
%Reference file (to align and flip reads)
codA = upper('catgaggggatgggtgcccgggtgaccgcatctcatacgactgccatgcattcttacaatggagcgtatacatcaagactgttccgcctgttgaagatgtccggtattaacttcgtggcgaaccccttggtcaatattcatttgcaggggagatttgatacttacccaaaaaggcgcgggataaccagggtaaaggagatgttggaatctggaattaacgtctgtttcggacacgacgacgtctttgacccc');
%Alignment threshold
thresh = 0.75;

%ORF reference files. EGFP, iRFP and mOrange should be present in the
%library, while BFP and mCherry should not

mCh = upper('gcaggacggcgagttcatctacaaggtgaagctgcgcggcaccaacttcccctccgacggccccgtaatgcagaagaaaaccatgggctggcaggcctcc');
EGFP = upper('tgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaacta');
iRFP = upper('TCCGGCATGGTCATCGGCGAAGCGAAGAGATCCGATCTGGAGTCATTCCTCGGCCAGCATTTCCCGGCATCTTTGGTGCCTCAGCAAGCGAGGTTGCTTT');
BFP = upper('aggccttcaccgagacactgtacccagccgatggtggactggagggcaggaatgacatggccctcaagctcgtcggaggcagtcacttgattgccaacgc');
mOrange = upper('CTATTTTAAGTTGAGCTTCCCAGAAGGATTTAAGTGGGAACGGGTGATGAATTACGAAGATGGGGGTGTAGTGACCGTTACTCAAGATAGCTCTTTGCAG');


reads = {};

%The total sequencing output from this experiment was > 20Gb, and will use
%a lot of memory if we load all the reads into one cell array (with
%fastqall). Instead, we'll analyze 5 fastq files at a time (in a for loop), align them to
%the reference and ORFs, split the reads into promoter and terminator
%regions (for all three ORFs), and discard reads from the original fastq
%file

for z = 1:5:length(dir)
    disp(strcat('Processing file number', num2str(z)))
    seq = cell(20000, 1);

%Load the files (5 at a time), and store the reads in an array

    for k = 1:5
        if direc(z + k - 1).isdir == 0
           n = direc(z + k - 1).name;
           [~, s, ~] = fastqread(n);
           seq((4000*k + 1):(4000*k+length(s)), 1) = s;
        end
    end

%filter reads by length
q = seq(~cellfun(@isempty, seq));

%Align to reference, flip as necessary
[~,~,pKR230,~] = alignup(q,codA,thresh);
%Remove reads that did not align to the reference
pKR230 = pKR230(~contains(pKR230, 'X'));

%Initialize vectors to store alignment percentages to the 5 ORFs
perc_mCh = zeros(length(pKR230), 1);
perc_gfp = zeros(length(pKR230), 1);
perc_iR = zeros(length(pKR230), 1);
perc_BFP = zeros(length(pKR230), 1);
perc_mOr = zeros(length(pKR230), 1);

start_points = zeros(length(pKR230), 5); %To store the start points

%Next, we'll align the reads to the 5 ORFs (3 of which should be there, and
%2 are negative controls that shouldn't be found)

for i = 1:length(pKR230)
%    if mod(i, 1000) == 0
%    disp(strcat('Aligning_', num2str(i)))
%    end
    
    a = cell2mat(pKR230(i));

%    [~, x1, y1] = swalign(a, mCh, 'alphabet', 'NT');
%    perc_mCh(i) = sum(x1(2, :) == '|')/length(mCh);
%    start_points(i, 1) = y1(1);
    [~, x1, y1] = swalign(a, EGFP, 'alphabet', 'NT');
    start_points(i, 2) = y1(1);
    perc_gfp(i) = sum(x1(2, :) == '|')/length(EGFP);
    [~, x1, y1] = swalign(a, iRFP, 'alphabet', 'NT');
    start_points(i, 3) = y1(1);
    perc_iR(i) = sum(x1(2, :) == '|')/length(iRFP);
%    [~, x1, y1] = swalign(a, BFP, 'alphabet', 'NT');
%    start_points(i, 4) = y1(1);
%    perc_BFP(i) = sum(x1(2, :) == '|')/length(BFP);
    [~, x1, y1] = swalign(a, mOrange, 'alphabet', 'NT');
    start_points(i, 5) = y1(1);
    perc_mOr(i) = sum(x1(2, :) == '|')/length(mOrange);

end

%Logical n X 5 array indicating whether each of the ORFs were found in each
%read (for n reads)

ORFs = [perc_mCh > 0.78, perc_gfp > 0.78, perc_iR > 0.78, perc_BFP > 0.78, perc_mOr > 0.78];

%Confusion matrix for the ORFs
%oconf = zeros(5, 5);
%for i = 1:5
%    for j = 1:5
%        oconf(i, j) = sum(ORFs(:, i) & ORFs(:, j));
%    end
%end

%figure
%confusionchart(oconf, {'mCherry', 'EGFP', 'iRFP670', 'BFP', 'mOrange'})

%Start points for the corresponding ORFs
start_points_eGFP = start_points(ORFs(:, 2) & ORFs(:, 3) & ORFs(:, 5), 2);
start_points_iRFP = start_points(ORFs(:, 2) & ORFs(:, 3) & ORFs(:, 5), 5);
start_points_mOr = start_points(ORFs(:, 2) & ORFs(:, 3) & ORFs(:, 5), 3);

%figure
%histogram(start_points_eGFP, 100, 'EdgeColor', 'None')
%hold on
%histogram(start_points_iRFP, 100, 'EdgeColor', 'None')
%hold on
%histogram(start_points_mOr, 100, 'EdgeColor', 'None')
%hold off


%Split Reads into promoter and terminator regions (for pKR230)
pKR230 = pKR230(ORFs(:, 2) & ORFs(:, 3) & ORFs(:, 5)); %Only keep reads that had all three ORFs

%Initialize vector to store the truncated regions of every read (promoter
%and terminator regions around all three ORFs

preads_eGFP = cell(size(start_points_eGFP));
preads_iRFP = cell(size(start_points_iRFP));
preads_mOr = cell(size(start_points_mOr));
treads_eGFP = cell(size(start_points_eGFP));
treads_iRFP = cell(size(start_points_iRFP));
treads_mOr = cell(size(start_points_mOr));


for i = 1:length(pKR230)
%    if mod(i, 100) == 0 %Uncomment to see progress (this part runs fairly quickly so shouldn't be necessary)
%        disp(i)
%    end
    a = cell2mat(pKR230(i)); %Load in the read
    
%If all three ORFs are present in the regions that we roughly expect them to be, split the read into the promoter and terminator regions    
    if start_points_eGFP(i) < 4000 & start_points_iRFP(i) < 6000 & start_points_mOr(i) < 8000
        preads_eGFP(i) = cellstr(a(1:start_points_eGFP(i)));
        treads_eGFP(i) = cellstr(a((start_points_eGFP(i)):start_points_iRFP(i)));
        preads_iRFP(i) = cellstr(a(start_points_eGFP(i):start_points_iRFP(i)));
        treads_iRFP(i) = cellstr(a((start_points_iRFP(i)):start_points_mOr(i)));
        preads_mOr(i) = cellstr(a(start_points_iRFP(i):start_points_mOr(i)));
        if length(a(start_points_mOr(i):end)) > 2000
            treads_mOr(i) = cellstr(a((start_points_mOr(i)):start_points_mOr(i)+2000));
        else
            treads_mOr(i) = cellstr(a(start_points_mOr(i):end));
        end
        
    else
        preads_eGFP(i) = cellstr('X');
        treads_eGFP(i) = cellstr('X');
        preads_iRFP(i) = cellstr('X');
        treads_iRFP(i) = cellstr('X');
        preads_mOr(i) = cellstr('X');
        treads_mOr(i) = cellstr('X');
    end    
end

%Append the 6 regions to the cell array "reads"
reads(end+ 1:end + length(preads_eGFP), 1:6) = [preads_eGFP, preads_mOr, preads_iRFP, treads_eGFP, treads_mOr, treads_iRFP];

end

%Split and store the promoter and terminator regions
preads_eGFP = reads(:, 1);
preads_mOr = reads(:, 2);
preads_iRFP = reads(:, 3);

treads_eGFP = reads(:, 4);
treads_mOr = reads(:, 5);
treads_iRFP = reads(:, 6);

%% Filter out reads that didn't contain all three expression units
x = contains(preads_eGFP, 'X') | contains(preads_iRFP, 'X') | contains(preads_mOr, 'X');

preads_eGFP = preads_eGFP(~x);
preads_iRFP = preads_iRFP(~x);
preads_mOr = preads_mOr(~x);

treads_eGFP = treads_eGFP(~x);
treads_iRFP = treads_iRFP(~x);
treads_mOr = treads_mOr(~x);

%Read promoter and terminator reference files

promoters = readcell('27k-Promoters.xlsx');
promoters = promoters(2:6, :);
promoters = promoters(:, 2);

terminators = readcell('27k-Terminators.xlsx');
terminators = terminators(:, 2);

%% Promoter and terminator assignments

%Promoter Assignments for pKReGFP (containment)

%Number of tiles contained for each promoter for every read (n X m matrix,
%where n is the number of reads, and m is the number of reference files)

ptiles_eGFP = zeros(length(preads_eGFP), 5);
ptiles_mOr = zeros(length(preads_mOr), 5);
ptiles_iRFP = zeros(length(preads_iRFP), 5);

for p = 1:5 %For every reference promoter
    disp(p)
    a = upper(cell2mat(promoters(p)));
    for m = 1:length(a) - 12
        b = a(m:m+11); %Break into 12bp tiles
        ptiles_eGFP(:, p) = ptiles_eGFP(:, p) + contains(preads_eGFP, b); %Add 1 for the reads that contained the tile
        ptiles_mOr(:, p) = ptiles_mOr(:, p) + contains(preads_mOr, b);
        ptiles_iRFP(:, p) = ptiles_iRFP(:, p) + contains(preads_iRFP, b);
    end
end

%Scale the number of tiles found by the length of the promoters, so every
%read is normalized against it (and represents what fraction of the tiles
%were found for each promoter)

promoter_lengths = repmat([length(cell2mat(promoters(1))), length(cell2mat(promoters(2))), length(cell2mat(promoters(3))), length(cell2mat(promoters(4))), length(cell2mat(promoters(5)))], length(preads_eGFP), 1);
ptiles_eGFP_scaled = ptiles_eGFP./promoter_lengths;
ptiles_mOr_scaled = ptiles_mOr./promoter_lengths;
ptiles_iRFP_scaled = ptiles_iRFP./promoter_lengths;

%Make EU specific confusion matrices for the promoters

pconf_eGFP = zeros(5, 5);
pconf_mOr = zeros(5, 5);
pconf_iRFP = zeros(5, 5);

for i = 1:5
    for j = 1:5
        %If a read contained at least 4% of all the tiles for a promoter,
        %we assign that promoter to the read
        pconf_eGFP(i, j) = sum(ptiles_eGFP_scaled(:, i) > 0.04 & ptiles_eGFP_scaled(:, j) > 0.04);
        pconf_mOr(i, j) = sum(ptiles_mOr_scaled(:, i) > 0.04 & ptiles_mOr_scaled(:, j) > 0.04);
        pconf_iRFP(i, j) = sum(ptiles_iRFP_scaled(:, i) > 0.04 & ptiles_iRFP_scaled(:, j) > 0.04);
    end
end

figure
subplot(1, 3, 1); confusionchart(pconf_eGFP); title('eGFP Promoters')
subplot(1, 3, 2); confusionchart(pconf_iRFP); title('iRFP Promoters')
subplot(1, 3, 3); confusionchart(pconf_mOr); title('mOrange Promoters')

%Containment for assigning terminators

ttiles_eGFP = zeros(length(treads_eGFP), 6);
ttiles_mOr = zeros(length(treads_mOr), 6);
ttiles_iRFP = zeros(length(treads_iRFP), 6);

for p = 1:6
    disp(p)
    a = upper(cell2mat(terminators(p)));
    for m = 1:length(a) - 11
        b = a(m:m+11);
        ttiles_eGFP(:, p) = ttiles_eGFP(:, p) + contains(treads_eGFP, b);
        ttiles_mOr(:, p) = ttiles_mOr(:, p) + contains(treads_mOr, b);
        ttiles_iRFP(:, p) = ttiles_iRFP(:, p) + contains(treads_iRFP, b);
    end
end

terminator_lengths = repmat([length(cell2mat(terminators(1))), length(cell2mat(terminators(2))), length(cell2mat(terminators(3))), length(cell2mat(terminators(4))), length(cell2mat(terminators(5))), length(cell2mat(terminators(6)))], length(treads_eGFP), 1);
ttiles_eGFP_scaled = ttiles_eGFP./terminator_lengths;
ttiles_mOr_scaled = ttiles_mOr./terminator_lengths;
ttiles_iRFP_scaled = ttiles_iRFP./terminator_lengths;

tconf_eGFP = zeros(6, 6);
tconf_mOr = zeros(6, 6);
tconf_iRFP = zeros(6, 6);

for i = 1:6
    for j = 1:6
        tconf_eGFP(i, j) = sum(ttiles_eGFP_scaled(:, i) > 0.03 & ttiles_eGFP_scaled(:, j) > 0.03);
        tconf_mOr(i, j) = sum(ttiles_mOr_scaled(:, i) > 0.03 & ttiles_mOr_scaled(:, j) > 0.03);
        tconf_iRFP(i, j) = sum(ttiles_iRFP_scaled(:, i) > 0.03 & ttiles_iRFP_scaled(:, j) > 0.03);
    end
end

%Terminators 1 and 2 are very similar, and thus have a lot of confusion
%between them. In the following section, we'll filter out reads that have
%>3% of tiles for either of these two terminators. If the number of tiles
%contained for any of the two terminators is at least 3% greater than the
%other, we assign that read to the corresponding terminator. Else we assign
%them to both/neither (both for the purposes of the confusion matrix, since
%we are unable to get a good assignment on those and need to discard them
%for further analysis downstream

x = ttiles_eGFP_scaled(ttiles_eGFP_scaled(:, 1) > 0.03 | ttiles_eGFP_scaled(:, 2) > 0.03, 1:2);
r1 = 0; r2 = 0; r3 = 0;
for k = 1:size(x, 1)
    if x(k, 1) > x(k, 2) + 0.03
        r1 = r1 + 1;
    elseif x(k, 2) > x(k, 1) + 0.03
        r2 = r2 + 1;
    else
        r3 = r3 + 1;
    end
end

tconf_eGFP(1, 1) = r1; tconf_eGFP(2, 2) = r2; tconf_eGFP(1, 2) = r3; tconf_eGFP(2, 1) = r3;

%Repeat for the other two ORFs

x = ttiles_mOr_scaled(ttiles_mOr_scaled(:, 1) > 0.03 | ttiles_mOr_scaled(:, 2) > 0.03, 1:2);
r1 = 0; r2 = 0; r3 = 0;
for k = 1:size(x, 1)
    if x(k, 1) > x(k, 2) + 0.03
        r1 = r1 + 1;
    elseif x(k, 2) > x(k, 1) + 0.03
        r2 = r2 + 1;
    else
        r3 = r3 + 1;
    end
end

tconf_mOr(1, 1) = r1; tconf_mOr(2, 2) = r2; tconf_mOr(1, 2) = r3; tconf_mOr(2, 1) = r3;

x = ttiles_iRFP_scaled(ttiles_iRFP_scaled(:, 1) > 0.03 | ttiles_iRFP_scaled(:, 2) > 0.03, 1:2);
r1 = 0; r2 = 0; r3 = 0;
for k = 1:size(x, 1)
    if x(k, 1) > x(k, 2) + 0.03
        r1 = r1 + 1;
    elseif x(k, 2) > x(k, 1) + 0.03
        r2 = r2 + 1;
    else
        r3 = r3 + 1;
    end
end

tconf_iRFP(1, 1) = r1; tconf_iRFP(2, 2) = r2; tconf_iRFP(1, 2) = r3; tconf_iRFP(2, 1) = r3;

%Plot confusion matrices
figure
subplot(1, 3, 1); confusionchart(tconf_eGFP); title('eGFP Terminators')
subplot(1, 3, 2); confusionchart(tconf_iRFP); title('iRFP Terminators')
subplot(1, 3, 3); confusionchart(tconf_mOr); title('mOrange Terminators')

%Since there seems to be almost no confusion in the promoter or terminator
%assignments, we can assign the promoter with the highest percentage of
%tiles as the promoter for that read. If a read had no tiles for any of the
%promoters or terminators, the "max" function will assign the first entry
%as the max by default. So let's add a column of all 0s at the start of the
%matrix, so that any assignments to column 1 would mean that the read had
%no promoter/terminator

%Assign promoter and terminator to an EU: iRFP
ptiles_eGFP_scaled(:, 6) = 0;
ptiles_eGFP_scaled(:, 2:6) = ptiles_eGFP_scaled(:, 1:5); %Move columns 1-5 to 2-6
ptiles_eGFP_scaled(:, 1) = 0; %Add all zeros in column 1

[~, prom_eGFP] = max(ptiles_eGFP_scaled');

ttiles_eGFP_scaled(:, 7) = 0;
ttiles_eGFP_scaled(:, 2:7) = ttiles_eGFP_scaled(:, 1:6);
ttiles_eGFP_scaled(:, 1) = 0;

[~, term_eGFP] = max(ttiles_eGFP_scaled');

%Number of promoter - terminator combinations for each ORF

assignment_table_eGFP = zeros(6, 7);
for i = 1:6
    for j = 1:7
        assignment_table_eGFP(i, j) = sum(prom_eGFP == i & term_eGFP == j);
    end
end

%Assign promoter and terminator to an EU: iRFP
ptiles_iRFP_scaled(:, 6) = 0;
ptiles_iRFP_scaled(:, 2:6) = ptiles_iRFP_scaled(:, 1:5);
ptiles_iRFP_scaled(:, 1) = 0;

[~, prom_iRFP] = max(ptiles_iRFP_scaled');

ttiles_iRFP_scaled(:, 7) = 0;
ttiles_iRFP_scaled(:, 2:7) = ttiles_iRFP_scaled(:, 1:6);
ttiles_iRFP_scaled(:, 1) = 0;

[~, term_iRFP] = max(ttiles_iRFP_scaled');

assignment_table_iRFP = zeros(6, 7);

for i = 1:6
    for j = 1:7
        assignment_table_iRFP(i, j) = sum(prom_iRFP == i & term_iRFP == j);
    end
end

%Assign promoter and terminator to an EU: mOr
ptiles_mOr_scaled(:, 6) = 0;
ptiles_mOr_scaled(:, 2:6) = ptiles_mOr_scaled(:, 1:5);
ptiles_mOr_scaled(:, 1) = 0;

[~, prom_mOr] = max(ptiles_mOr_scaled');

ttiles_mOr_scaled(:, 7) = 0;
ttiles_mOr_scaled(:, 2:7) = ttiles_mOr_scaled(:, 1:6);
ttiles_mOr_scaled(:, 1) = 0;

[~, term_mOr] = max(ttiles_mOr_scaled');

assignment_table_mOr = zeros(6, 7);

for i = 1:6
    for j = 1:7
        assignment_table_mOr(i, j) = sum(prom_mOr == i & term_mOr == j);
    end
end

%% Downstream analysis (post promoter/terminator assignment)

%Assigning an EU ARRAY (promoter terminator combinations for each ORF) for
%every read

EU_array = zeros(size(prom_eGFP));

for r = 1:length(preads_eGFP)
    if prom_eGFP(r) > 1 && prom_mOr(r) > 1 && prom_iRFP(r) > 1 && term_eGFP(r) > 1 && term_mOr(r) > 1 && term_iRFP(r) > 1
        EU_array(r) = (prom_eGFP(r)-2)*5400 + (term_eGFP(r)-2)*900 + (prom_iRFP(r)-2)*180 + (term_iRFP(r)-2)*30 + (prom_mOr(r)-2)*6 + term_mOr(r) - 1;
    end
end

line = 1:1:27000;
array_numbers = zeros(size(line));

for i = 1:27000
    disp(i)
    array_numbers(i) = sum(EU_array == i);
%    array_numbers(i) = 0;
end

%Find number of EUs found at lower read depths

read_depth = 1:1000:length(EU_array);
eus_found = zeros(size(read_depth));

for k = 1:length(read_depth)
    y = read_depth(k);
%    if mod(k, 1000) == 0
    disp(k)
%    end
    x = datasample(EU_array, y, 'Replace',false);
    eus_found(k) = length(unique(x)) - 1;
end







