clear all;
%%TASK 1 investigate animal:Artic fox(common name):Vulpes lagopus(taxonmic name) 
%%Accession number:NC_026529
%%cytochrome B accession number:YP_009122397
%%task 1 :statistical analysis
close all;
seq = getgenbank('NC_026529','SequenceOnly',true);
length(seq)
basecount(seq)
ntdensity(seq,'window',2000)
ntdensity(seq,'window',500)
[dimers,matrix]=dimercount(seq)

%set the minimum orflength, and find the ORFS
artic_orfs=seqshoworfs(seq,'geneticcode',2,'frames','all','minimumlength',379,'nodisplay','on')
o=0;
for k=1:6
    if (length(artic_orfs(k).Stop))==length(artic_orfs(k).Start)
        o=[o,artic_orfs(k).Stop(1:end)-artic_orfs(k).Start(1:end)];
    else
        o=[o,artic_orfs(k).Stop(1:end)-artic_orfs(k).Start(1:end-1)]
    end
end
o=o(2:end)

%%find protein-coding genes, and translate CYTB and CYTC
seq_cytB = getgenpept('YP_009122397','SequenceOnly',true);
seq_cytC = getgenpept('YP_009122391','SequenceOnly',true);
orfs=seqshoworfs(seq,'geneticcode',2,'frames','all','minimumlength',64,'nodisplay','on');
index =0;
for k=1:6
    a= orfs(k).Stop;
    b= orfs(k).Start;
    if length(a)==length(b)
        c = (a-b)/3
         index = find(c==379)
         if index~=0
    
             frameValue = k;
             index = index;
         end
    else
        continue;
    end
end
cytb_orf = seq(orfs(2).Start(4):orfs(2).Stop(4))
%%translate nucleotides to AA
cytb_seq = nt2aa(cytb_orf,'GeneticCode',2)
%%to find cytoC
for k=1:6
    a= orfs(k).Stop;
    b= orfs(k).Start;
    if length(a)==length(b)
        c = (a-b)/3
         index = find(c==227)
         if index~=0
             frameValue = k;
             index = index;
         end
    else
        continue;
    end
end
cytc_orf = seq(orfs(1).Start(1):orfs(1).Stop(1))
%%translate nucleotides to AA
cytc_seq = nt2aa(cytb_orf,'GeneticCode',2)
sc1 = seqshoworfs(cytc_seq);
sc2 = seqshoworfs(seq_cytC);
% Task 3/ using grey wolf as outgroup
seqs=fastaread('coursework.fasta')
dist=seqpdist(seqs,'method','jukes-cantor','indels','pair','Alphabet','AA');
tree=seqneighjoin(dist,'equivar',seqs)
view(tree)
newtree=reroot(tree,9)
view(newtree)

%Task 4 AND 5 add gray wolf, Golden jackal and coyote, using polar bear as
%outgroup
%% cytochrome B/C AA
seqs_4b=fastaread('add_new_species (1).fasta');
distb=seqpdist(seqs_4b,'method','jukes-cantor','indels','pair');
tree_BAA=seqneighjoin(distb,'equivar',seqs_4b)
view(tree_BAA)
seqs_4c=fastaread('cytoc_task4.fasta');
distc=seqpdist(seqs_4c,'method','jukes-cantor','indels','pair');
tree_CAA=seqneighjoin(distc,'equivar',seqs_4c)
newtree_aac=reroot(tree_CAA,4)
plot(newtree_aac)
%% Cyt B NT
seqsbnt = fastaread('cytoBNT.fasta');
distb_nt = seqpdist(seqsbnt,'method','jukes-cantor','indels','pair','alphabet','NT');
tree = seqneighjoin(distb_nt,'equivar',seqsbnt);
newtree_bnt = reroot(tree,5);
view(newtree_bnt)
%% Cyt C NT
seqs = fastaread('cytoCNT.fasta');
dist_cnt = seqpdist(seqs,'method','jukes-cantor','indels','pair','alphabet','NT');
tree = seqneighjoin(dist_cnt,'equivar',seqs);
newtree_cnt = reroot(tree,3);
view(newtree_cnt)
%consensus tree
allNames = {'artic fox','kit fox','red fox','golden jackal','gray wolf','coyote','polar bear'};
weights = [sum(distb) sum(distc) sum(distb_nt) sum(dist_cnt)];
weights = weights / sum(weights);
dist_con = distb .*weights(1)+distc .*weights(2)+distb_nt .*weights(3)+dist_cnt .*weights(4);
tree_con = seqlinkage(dist_con,'average',allNames);
plot(tree_con,'type','angular');
title('consensus tree')
%%task 7 alignment
% alignment cytochrome B
ali1 = fastaread('add_new_species (1).fasta')
ma = multialign(ali1,'verbose',true)
seqalignviewer(ma)
%alignment cytochrome C
ali2 = fastaread('cytoc_task4.fasta')
ma2 = multialign(ali2,'verbose',true)
seqalignviewer(ma2)
