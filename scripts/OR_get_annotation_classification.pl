#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# Perl script to categorize the gene annotations as complete, partial/pseudogenes or fragment based in the the homology-based alignments (BLASTP and HMMER), and the metrics set by the user

my $bitacorapath = "/sietch_colab/scamison/gen-families/bitacora-BITACORA-v1.4.2"; # Specify the path to bitacora folder

# Input files | Specify the path and file for each of the required input files

my $genome = "data/culexPipiens.softmasked.fasta "; # Genome assembly in FASTA
my $annotgff = "OR_PF02949/OR_PF02949_genomic_and_annotated_genes.gff3"; # Gene family annotation GFF3 file
my $annotprot = "OR_PF02949/OR_PF02949_genomic_and_annotated_proteins_trimmed.fasta"; # Gene family annotation protein sequences in FASTA file

my $genefamilyname = "OR"; # Short name of the gene family to use in the output files - Do not use spaces here
my $blastoutput = "OR_PF02949/OR_PF02949_genomic_and_annotated_proteins_trimmed.fasta_ORblast.txt"; # BLAST output aligning the annotated gene family
my $filtergenefamilyname = "GR"; # Short name of the gene family or families that are used to exlude incorrect annotations (i.e. ORs from GR annotations)
my $filterblastoutput = "OR_PF02949/OR_PF02949_genomic_and_annotated_proteins_trimmed.fasta_GRblast.txt"; # BLAST output aligning the gene family to exclude (i.e. ORs from GR annotations)
my $hmmeroutput = "OR_PF02949/OR_PF02949_genomic_and_annotated_proteins_trimmed.fasta_hmmer.domtblout"; # HMMER output aligning the protein domain into the gene family annotations

my $dmelblastoutput = "OR_PF02949/OR_PF02949_genomic_and_annotated_proteins_trimmed.fasta_ORcleanblast.txt"; # BLAST output aligning the annotated gene family from Drosophila melanogaster


# Length parameters | Specify here the required length of the homology-based alignments to classify the genes into complete, partial or fragments

my $completelength = "302"; # Length to consider complete gene copies
my $minlength = "175"; # Minimum length to retain partial/pseudogenes - Seqs shorter than "minlength" but longer than 50aa will be labelled as fragment
my $maxlength = "477"; # Maximum length, it will print a warning with the sequences larger than this length, which could be chimeric genes.


# Code

my ($line, $name, $nameout);

my $prefix = "";
if ($annotgff =~ /(\S+)\.gff3/){
	$prefix = $1;
} else {die "Can't find prefix on gff3 file $annotgff\n";}

my $outgffall = "$prefix"."\_OR_all.gff3"; 
my $outprotall = "$prefix"."\_OR_all"; 

my $outgff = "$prefix"."\_OR_all_nofragment.gff3"; 
my $outprot = "$prefix"."\_OR_all_nofragment"; 

my $outgffcomp = "$prefix"."\_OR_complete.gff3"; 
my $outprotcomp = "$prefix"."\_OR_complete"; 


# Reading Protein fasta

my %fasta;
open(File, "<", $annotprot) or die "Cannot open file '$annotprot': $!\n";
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	
	$line =~ s/_\d+dom//; # Remove the _Ndom tag from the genes, which is not present in the gff3 file
	if ($line =~ />(\S+)/){
		$name = $1;
	} else {
		$fasta{$name} .= "$line";
	}
}
close File;


# Reading hmmer output

my $evalue = "0.001";
my %grhmmer; my %hmmer;
open (File, "<", $hmmeroutput) or die "Cannot open file '$hmmeroutput': $!\n";
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
    next if ($line =~ /^#/);
    $line =~ s/\s+/\t/g;$line =~ s/\t+/\t/g;
    $line =~ s/_\d+dom//; # Remove the _Ndom tag from the genes, which is not present in the gff3 file
    my @subline = split (/\t/, $line);
    if ($subline[6] <= $evalue){ # E-value
            next if ($subline[12] > $evalue); # E-value specific domain region
            push (@{$hmmer{$subline[0]}}, join("\t",$subline[0], $subline[17], $subline[18]));
            #$multipledom{$subline[0]}++;
            #$length{$subline[0]} = $subline[2];
            $grhmmer{$subline[0]}++;
    }
}
close File;


# Reading blast output

my %orblast;
open (File, "<", $filterblastoutput) or die "Cannot open file '$filterblastoutput': $!\n";
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	
	$line =~ s/_\d+dom//; # Remove the _Ndom tag from the genes, which is not present in the gff3 file
	my @subline = split (/\t/, $line);

	if ($subline[10] <= $evalue){
		if (exists $orblast{$subline[0]}){
			if ($subline[10] < $orblast{$subline[0]}){
				$orblast{$subline[0]} = "$subline[10]";
			}
		} else {
			$orblast{$subline[0]} = "$subline[10]";
		}

	}
}
close File;


my %dmelblast; my %dmelblasteval;
open (File, "<", $dmelblastoutput) or die "Cannot open file '$dmelblastoutput': $!\n";
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	
	$line =~ s/_\d+dom//; # Remove the _Ndom tag from the genes, which is not present in the gff3 file
	my @subline = split (/\t/, $line);

	if ($subline[10] <= "1e-40"){
		if (exists $dmelblast{$subline[0]}){
			if ($subline[10] < $dmelblasteval{$subline[0]}){
				$dmelblasteval{$subline[0]} = "$subline[10]";
				$dmelblast{$subline[0]} = "$subline[1]";
			}
		} else {
			$dmelblasteval{$subline[0]} = "$subline[10]";
			$dmelblast{$subline[0]} = "$subline[1]";
		}

	}
}
close File;


my %blast; my %length;
my %grblast;
open (File, "<", $blastoutput) or die "Cannot open file '$blastoutput': $!\n";
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	
	$line =~ s/_\d+dom//; # Remove the _Ndom tag from the genes, which is not present in the gff3 file
	my @subline = split (/\t/, $line);

	if ($subline[10] <= $evalue){
		if (exists $grblast{$subline[0]}){
			if ($subline[10] < $grblast{$subline[0]}){
				$grblast{$subline[0]} = "$subline[10]";
			}
		} else {
			$grblast{$subline[0]} = "$subline[10]";
		}

		push (@{$blast{$subline[0]}}, join("\t",$subline[1], $subline[10], $subline[3], $subline[12], $subline[13], $subline[6], $subline[7], $subline[2])); # New line to save the whole blast results and explore for chimeric genes
		$length{$subline[0]} = $subline[12];
	}
}
close File;


### Finding chimeric/fused genes

my $minlengthcut = "30"; ## Minimum positions required to trim a protein (i.e. blast hits starting in position 10 will report the full sequence instead of trimming the first 10 aa)

foreach my $key (sort keys %blast) {
	my $blasthit2="";
	my $hitlvl = "0";
	my (@ini, @fin);
	$ini[0] = "99999999999999999999";
	$fin[0] = "1";
	foreach my $blastresult (@{$blast{$key}}) {
		my @subline = ();
		@subline = split (/\t/, $blastresult);
		my $hitblast = $key;
		my $filtro1 = ($subline[4]*2)/3;
		my $filtro2 = ($subline[3]*0.8);

#### NO FILTERING HERE		
#		if ($subline[2] < 51) { # If the alignment is lower than 50aa, and smaller than 2/3 QUERY length protein used, it should contain a similarity higher than 80%. If not it is removed as false positive (small domains hitting non-related proteins)
#			unless ($subline[2] >= $filtro1){
#				next if ($subline[7] < 80);
#			}
#		}

#		if ($subline[2] >= $filtro1 || $subline[2] >= $filtro2) { # BLAST filtering: alignment covering 2/3 of subject, or 80% of query
			$hitlvl = "2";
			my $n = 0;
			my $extrahit = 0;
			foreach my $i (@ini){
				my $f = $fin[$n];
				my $f2 = $f + 10;
				my $i2 = $i - 10;

				if ($i >= int($subline[5]) && $f <= int($subline[6])){
					$ini[$n] = int($subline[5]);
					$fin[$n] = int($subline[6]);
				}
				elsif ($i <= int($subline[5]) && $f < int($subline[6]) && $f2 >= int($subline[5])){
					$fin[$n] = int($subline[6]);
				}
				elsif ($i > int($subline[5]) && $f >= int($subline[6]) && $i2 <= int($subline[6])){
					$ini[$n] = int($subline[5]);
				}
				elsif ($f2 < int($subline[5]) || $i2 > int($subline[6])) {
					$extrahit++;
				}

				$n++;
			}

			if ($extrahit >= $n) {
					$ini[$n] = int($subline[5]);
					$fin[$n] = int($subline[6]);				
			}
#		}
	}

	if ($hitlvl == 2){
		my $hits = scalar(@ini);

		if ($hits == 1){

			# Length filter to avoid exluding a few ($minlengthcut) initial or end positions
			my $ipos = $ini[0];
			my $fpos = $fin[0];
			if ($ipos <= $minlengthcut){ # Initial position
				$ipos = 1;
			}
			my $filterend = $length{$key} - $minlengthcut;
			if ($fpos >= $filterend){ # Initial position
				$fpos = $length{$key};
			}

#			print Results "$key annot $ipos $fpos blastp\n";
		}
		else {
			my $n = 0;

			print "Warning $key in with length $length{$key} should be split in $hits parts\n";
			
			foreach my $i (@ini){
				my $f = $fin[$n];

				# Length filter to avoid exluding a few ($minlengthcut) initial or end positions
				my $ipos = $ini[$n];
				my $fpos = $fin[$n];
				if ($ipos <= $minlengthcut){ # Initial position
					$ipos = 1;
				}
				my $filterend = $length{$key} - $minlengthcut;
				if ($fpos >= $filterend){ # Initial position
					$fpos = $length{$key};
				}

				my $nn= $n+1;
#				print Results "$key\_split$nn annot $ipos $fpos blastp\n";
				$n++;

			}
		}
	}

}

###


# Reading GFF to obtain gene order and exon number

my %clusterorder;
my %exonnumber;
open (GFFfile , "<", $annotgff) or die "Cannot open file '$annotgff': $!\n";
while (<GFFfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	next if ($line =~ /^#/);
	my @subline = split (/\t/, $line);


	if ($subline[2] =~ /mRNA/){

		my $genename = "";
		if ($subline[8] =~ /ID=([^;]+)/){
			$genename = $1;
		}
		else {print "ERROR in run_OR_classification: It fails detecting ID in GFF file $line\n";}

#		$clusterorder{$seannameclust}{$seannameclustcoord}{$genenamesean} = $genename; 
		$clusterorder{$subline[0]}{$subline[3]}{$genename} = $genename; 


	} elsif ($subline[2] =~ /CDS/){

		my $genename = "";
		if ($subline[8] =~ /Parent=([^;]+)/){
			$genename = $1;
		}
		else {print "ERROR in run_OR_classification: It fails detecting ID in GFF file $line\n";}

		$exonnumber{$genename}++;

	}
}
close GFFfile;


# Creating table and dividing proteins

my %newnames;
my $deleteproteins = "";
my $deleteproteinsl50 = "";
my %hstatus;

open (Results, ">", "$genefamilyname\_gene_table.tsv");
print Results "Gene name\tScaffold location\tFinal gene name\tLength\tExon number\tBlast hit\tDomain\tDmel function\tStatus Complete Partial Pseudogene\tFragment\n";

open (Resultssum, ">", "$genefamilyname\_summary_table.tsv");
#print Resultssum "$genefamilyname total number\tComplete\tComplete with domain\tPartial\tPartial with domain\tPartial longer 350aa\tPseudogene\tPseudogene with domain\tPseudogene longer 350aa\tLength lower than 200aa\tLength lower than 100aa\tLength lower than 50aa\tOR count\tOrco\tOrco complete\tSingle-exon GRs complete\tSingle-exon GRs incomplete\tSingle-exon GRs incomplete (length > 350aa)\tSingle-exon GRs pseudogenes\tDmel sugar receptors\tDmel bitter receptors\tDmel CO2 receptors\tDmel fructose receptors\tGR complete + partial with domain\tFinal retained GRs\tFinal retained complete GRs (complete + partial >350aa)\tFinal retained complete de novo annot\tFinal retained partial GRs\tFinal retained partial de novo annot\tFinal retained pseudogene GRs\tFinal retained complete 1-exon GRs\tFinal retained incomplete or pseudogene 1-exon GRs\tFragment GRs partial\tFragment GRs pseudogenes\n";
#print Resultssum "$totalor\t$com\t$comd\t$par\t$pard\t$parc\t$pse\t$psed\t$psec\t$l200\t$l100\t$l50\t$orcount\t$orconum\t$orconumcom\t$ninec\t$ninei\t$nineic\t$ninep\t$sugar\t$bitter\t$co2\t$gr43\t$compard\t$ftot\t$fcomp\t$denovoc\t$fpar\t$denovop\t$fpse\t$fninec\t$fninei\t$frpartial\t$frpse\n";
print Resultssum "$genefamilyname total annotated genes\t$filtergenefamilyname genes excluded\tGenes shorter than 50aa excluded\tFinal retained $genefamilyname\tComplete $genefamilyname\tComplete annotated de novo\tPartial $genefamilyname\tPseudogenes\tFragment copies\n";
#print Resultssum "$totalor\t$orcount\t$l50\t$ftot\t$fcomp\t$denovoc\t$fpar\t$fpse\t$frpartpse\n";

my $ornum = "1";
my $fragornum = "1";
my $totalor = "0"; my $com = "0"; my $comd = "0"; my $par = "0"; my $pard = "0"; my $parc = "0"; my $pse = "0"; my $psed = "0"; my $psec = "0"; my $l50 = "0"; my $l100 = "0"; my $l200 = "0"; my $orconum = "0"; my $orconumcom = "0"; my $ninec = "0"; my $ninei = "0"; my $nineic = "0"; my $ninep = "0";
my $compard = "0"; my $ftot = "0"; my $fcom = "0"; my $fpar = "0"; my $fparc = "0"; my $fcomp = "0"; my $fpse = "0"; my $fninec = "0"; my $fninei = "0"; my $frpartial = "0"; my $frpse = "0"; my $frpartpse = "0";
my $orcount = "0";
my $bitter = "0"; my $sugar = "0"; my $co2 = "0"; my $gr43 = "0";
my $denovoc = "0"; my $denovop = "0";

my %lsfasta;
my %lengthseq;
foreach my $prot (sort keys %fasta){
	my $length = length ($fasta{$prot});
	$lsfasta{$length}{$prot} = "$fasta{$prot}";
	$lengthseq{$prot} = $length;
}

#foreach my $len (sort { $b <=> $a } (keys %lsfasta)){
#	foreach my $prot (keys %{$lsfasta{$len}}){
foreach my $scaf (sort (keys %clusterorder)){
foreach my $coord (sort { $a <=> $b } (keys %{$clusterorder{$scaf}})){
	foreach my $genename (keys %{$clusterorder{$scaf}{$coord}}){

		my $prot = $clusterorder{$scaf}{$coord}{$genename};

		my $len = "0";
		my $seq = "";
		if (exists $lengthseq{$prot}){
			$len = $lengthseq{$prot};
			$seq = $lsfasta{$len}{$prot};
		} else { # Avoid errors in empty genes from genewise
#			die "Error: Can't find length of $prot, located in $scaf $coord $genename\n";
			$seq = "NNNNNNNNNNNNNXNNNXNNNNNNNNN";
			$exonnumber{$prot} = 0;
		}

		$totalor++;

		my $orcohit = "None";
=h		
		if (exists $orco{$prot}){
			if ($orco{$prot} =~ /Orco/){
				$orcohit = "Yes";
				$orconum++;
#				if ($status =~ /Complete/){
#					$orconumcom++;
#				}
			}
		}
=cut

		my $status = "";
		my $sufix = "";
		my $discard = "No"; # Label for fragment ORs

		if ($len >= $completelength){
			$status = "Complete";
			$hstatus{$prot} = "Complete";
#			$com++;
		} else {
			$status = "Partial";
			$hstatus{$prot} = "Partial";
			$sufix = "S";
#			$par++;
			if ($len < 200){
				$l200++;
#				$deleteproteins .= "$prot ";	
				#$discard = "Yes";				
			}
			if ($len < 100){
				$l100++;
#				$deleteproteins .= "$prot ";	
				#$discard = "Yes";				
			}
			if ($len < 50){
				$l50++;
#				$deleteproteins .= "$prot ";   # Delete all protein shorter than 50aa
				$deleteproteinsl50 .= "$prot ";   # Delete all protein shorter than 50aa
#				$discard = "Yes";
			}
			if ($len < $minlength){
				if ($orcohit =~ /None/){ # Avoid deleting ORCO even if it is partial
					$deleteproteins .= "$prot ";	# Delete all proteins shorter than the specified minimum length
					$discard = "Yes";		
				}		
			}				
		} 

		my $middleseq = substr ($seq, 5, -5);
		if ($middleseq =~ /X/){
			$status = "Pseudogene";
			$hstatus{$prot} = "Pseudogene";
			$sufix = "P";
#			$pse++;
		}


		my $orhit = "None";
		if (exists $orblast{$prot} && exists $grblast{$prot} ){
			if ($orblast{$prot} <= $grblast{$prot}){
				$orhit = "OR";
			} elsif ($orblast{$prot} > $grblast{$prot}){
				$orhit = "GR";
				$deleteproteinsl50 .= "$prot "; ### Delete protein identified as GR
				$discard = "Yes";
				$orcount++;
			}  
		} elsif (exists $orblast{$prot}){
			$orhit = "OR";
		} elsif (exists $grblast{$prot}){
			$orhit = "GR";
			$deleteproteinsl50 .= "$prot "; ### Delete protein identified as GR
			$discard = "Yes";
			$orcount++;
		}

		my $ordomain = "None";
		if (exists $grhmmer{$prot}){
			$ordomain = "OR";
			if ($status =~ /Pseudogene/){
				$psed++;
			} elsif ($status =~ /Complete/){
				$comd++;
			} elsif ($status =~ /Partial/){
				$pard++;
				$compard++;
			}
		}


		if ($orhit =~ /OR/){
#			$deleteproteinsl50 .= "$prot "; ### Delete protein identified as OR
#			$discard = "Yes";
#			$orcount++;
		} elsif ($ordomain =~ /OR/){
			if ($orhit !~ /GR/){
#				$deleteproteinsl50 .= "$prot "; ### Delete protein identified as OR
#				$discard = "Yes";
#				$orcount++;
			}
			# If blast hits OR, it is ok as some domain are mis-identified
		}


		if ($status =~ /Complete/){
			$com++;
			$compard++;
		} elsif ($status =~ /Partial/){
			$par++;
			if ($len >= $completelength){
				$status = "Complete"; ## Uncomment if we finally keep these also as complete copies
				$hstatus{$prot} = "Complete";
				$parc++;
			}

		} elsif ($status =~ /Pseudogene/){
			$pse++;
			if ($len >= $completelength){
				#$status = "Complete"; ## Uncomment if we finally keep these also as complete copies (Likely not as putative pseudogenes)
				#$hstatus{$prot} = "Complete";
				$psec++;
			}
		}		

		# Final filter to discard proteins
		if ($status =~ /Partial/ || $status =~ /Pseudogene/){
			if ($ordomain =~ /None/ && $len < $minlength){
				if ($orcohit =~ /None/){ # Avoid deleting ORCO even if it is partial and have no domain
					$deleteproteins .= "$prot "; ### Delete protein partial without domain, and length < specified in minorlength
					$discard = "Yes";		
				}	
			}
		}

		# Count 9 exon ORs || Single exon GRs 
		if ($exonnumber{$prot} == 1){
			if ($status =~ /Complete/){
				$ninec++;
			} elsif ($status =~ /Partial/) {
				$ninei++;
				if ($len >= $completelength){
					$nineic++;
				}				
			} elsif ($status =~ /Pseudogene/) {
				$ninep++;
			}
		}
		unless ($discard =~ /Yes/){
			if ($exonnumber{$prot} == 1){
				if ($status =~ /Complete/){
					$fninec++;
				} else {
					$fninei++;
				}	
			}
		}	


		# Find Dmel function 
		my $dmelhit = "None";
		if (exists $dmelblast{$prot}){
			if ($dmelblast{$prot} =~ /Gr64a/ || $dmelblast{$prot} =~ /Gr64b/ || $dmelblast{$prot} =~ /Gr64c/ || $dmelblast{$prot} =~ /Gr64d/ || $dmelblast{$prot} =~ /Gr64e/ || $dmelblast{$prot} =~ /Gr64f/ || $dmelblast{$prot} =~ /Gr5a/ || $dmelblast{$prot} =~ /Gr61a/){
				$dmelhit = "$dmelblast{$prot} Sugar";
				$sugar++;
			} elsif ($dmelblast{$prot} =~ /Gr21a/ || $dmelblast{$prot} =~ /Gr63a/){
				$dmelhit = "$dmelblast{$prot} CO2";
				$co2++;
			} elsif ($dmelblast{$prot} =~ /Gr43a/){
				$dmelhit = "$dmelblast{$prot} Fructose";
				$gr43++;
			} elsif ($dmelblast{$prot} =~ /Gr93a/ || $dmelblast{$prot} =~ /Gr66a/ || $dmelblast{$prot} =~ /Gr33a/){
				$dmelhit = "$dmelblast{$prot} Bitter";
				$bitter++;
			} else {
				$dmelhit = "$dmelblast{$prot}";
			}
		}


		# Get new names
		my $nname = "";
		$nname = $prot; # Skip here the renaming, use the gene names from the annotation | uncomment the following lines to rename gene ids

		if ($discard =~ /Yes/){
			if ($deleteproteinsl50 =~ /$prot /){
				$nname = "Discarded";
			}
		}

=head		
		if ($discard =~ /Yes/){
			#$nname = "Discarded";
			#if ($len < 50){
			#	$nname = "Discarded";
			if ($deleteproteinsl50 =~ /$prot /){
				$nname = "Discarded";
			} else {
				$nname = "$gagashort{$gagaid}"."fragGR"."$fragornum"."$sufix";
				$fragornum++;
			}
		} else {
			if ($orcohit =~ /Yes/){
				$nname = "$gagashort{$gagaid}"."Orco"."$sufix";
				if ($status =~ /Complete/){
					$orconumcom++;
				}
				if ($orconum > 1){
					$nname = "$gagashort{$gagaid}"."Orco"."$orconum"."$sufix";
				}
				$deleteproteinsl50 .= "$prot "; # Delete ORco from final retained GRs
			} else {
				$nname = "$gagashort{$gagaid}"."GR"."$ornum"."$sufix";
				$ornum++;
			}
			if ($len >= $maxlength){ # Print a warning in retained GRs longer than the max specified length
				print "Warning: $prot in $gagaid with length $len is longer than the max specified length\n";
			}
		}
=cut		
		$newnames{$prot} = "$nname";


		print Results "$prot\t$scaf\-$coord\t$nname\t$len\t$exonnumber{$prot}\t$orhit\t$ordomain\t$dmelhit\t$status\t$discard\n";

	}

}}
close Results;


# Reading GFF and generating final file filtering short proteins and renaming ORs

open (Resultsgff, ">", "$outgff");
open (Resultsgffcomp, ">", "$outgffcomp");
open (Resultsgffall, ">", "$outgffall");

open (GFFfile , "<", $annotgff); 
while (<GFFfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);

	if ($line =~ /^#/){
		print Resultsgff "$line\n";
		next;
	}
	
	my @subline = split (/\t/, $line);

	if ($subline[2] =~ /CDS/){
		my $genename = "";
		my $rest = "";
		if ($subline[8] =~ /Parent=([^;]+)(\S*)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
			$rest = $2;
		}
		else {die "ERROR in run_OR_classification.pl: It fails detecting Parent ID in $line\n";}

#		my $nnamef = $newnames{$genename};
		my $nnamef = "";
		if (exists $newnames{$genename}){
			$nnamef = $newnames{$genename};
		} else {
			die "Can't find gene new name for $genename in $line\n";
		}

# 		Keeping now all genes as fragment ORs; and creating multiple GFFs
#		next if ($deleteproteins =~ /$genename /); ### Delete genes that do not pass the filter
		next if ($deleteproteinsl50 =~ /$genename /); ### Delete from GFF all genes shorter than 50aa or mis-identification such as GR

		print Resultsgffall "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tParent=$nnamef"."$rest\;Status=$hstatus{$genename}\n";

		unless ($deleteproteins =~ /$genename /){
			print Resultsgff "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tParent=$nnamef"."$rest\;Status=$hstatus{$genename}\n";

			if ($hstatus{$genename} =~ /Complete/){
				print Resultsgffcomp "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tParent=$nnamef"."$rest\;Status=$hstatus{$genename}\n";
			}
		}
		
	}
	elsif ($subline[2] =~ /mRNA/){
		my $genename = "";
		my $rest = "";
		if ($subline[8] =~ /ID=([^;]+)\;Parent=[^;]+(\;\S*)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
			$rest = $2;
		}
		else {print "ERROR in run_OR_classification.pl: It fails detecting ID in $line\n";}

		my $nnamef = "";
		if (exists $newnames{$genename}){
			$nnamef = $newnames{$genename};
		} else {
			die "Can't find gene new name for $genename in $line\n";
		}


# 		Keeping now all genes as fragment ORs; and creating multiple GFFs
#		next if ($deleteproteins =~ /$genename /); ### Delete genes that do not pass the filter
		next if ($deleteproteinsl50 =~ /$genename /); ### Delete from GFF all genes shorter than 50aa or mis-identification such as GR

		print Resultsgffall "$subline[0]\t$subline[1]\tgene\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";		
		print Resultsgffall "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=$nnamef\;Parent=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";		

		unless ($deleteproteins =~ /$genename /){
			print Resultsgff "$subline[0]\t$subline[1]\tgene\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";		
			print Resultsgff "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=$nnamef\;Parent=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";		

			if ($hstatus{$genename} =~ /Complete/){
				print Resultsgffcomp "$subline[0]\t$subline[1]\tgene\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";		
				print Resultsgffcomp "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=$nnamef\;Parent=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";		
			}
		}

		
		# Obtain final OR number statistics
		if ($deleteproteins =~ /$genename /){
			$ftot++; # New, also counting fragments into the final retained counts, then specified how many fragments
			$frpartpse++;
			if ($hstatus{$genename} =~ /Complete/){
				# Nothing here
			} elsif ($hstatus{$genename} =~ /Partial/){
				$frpartial++;
			} elsif ($hstatus{$genename} =~ /Pseudogene/){
				$frpse++;
			}		

		} else {
			$ftot++;
			if ($hstatus{$genename} =~ /Complete/){
				$fcom++;
				$fcomp++;
				if ($genename =~ /GR.*\.t/){
					$denovoc++; # Counting new annotated genes in bitacora gemoma
				}
			} elsif ($hstatus{$genename} =~ /Partial/){
				$fpar++;
				if ($genename =~ /GR.*\.t/){
					$denovop++; # Counting new annotated genes in bitacora gemoma
				}
				if ($lengthseq{$genename} >= $completelength){ # It will not do nothing here as I relabelled those as complete in the status
					$fparc++;
					$fcomp++;
				}
			} elsif ($hstatus{$genename} =~ /Pseudogene/){
				$fpse++;
			}
		}

	}
}
close GFFfile;


#Encode proteins and CDS from the generated GFFs

system ("perl $bitacorapath\/Scripts/Tools/gff2fasta_v3.pl $genome $outgff $outprot ");
system ("sed \'s\/X\*\$\/\/\' $outprot\.pep.fasta > $outprot\.pep.fasta.tmp");
system ("mv $outprot\.pep.fasta.tmp $outprot\.pep.fasta");

system ("perl $bitacorapath\/Scripts/Tools/gff2fasta_v3.pl $genome $outgffall $outprotall ");
system ("sed \'s\/X\*\$\/\/\' $outprotall\.pep.fasta > $outprotall\.pep.fasta.tmp");
system ("mv $outprotall\.pep.fasta.tmp $outprotall\.pep.fasta");

system ("perl $bitacorapath\/Scripts/Tools/gff2fasta_v3.pl $genome $outgffcomp $outprotcomp ");
system ("sed \'s\/X\*\$\/\/\' $outprotcomp\.pep.fasta > $outprotcomp\.pep.fasta.tmp");
system ("mv $outprotcomp\.pep.fasta.tmp $outprotcomp\.pep.fasta");

#print Resultssum "$gagaid\t$gagasp{$gagaid}\t$totalor\t$com\t$comd\t$par\t$pard\t$parc\t$pse\t$psed\t$psec\t$l200\t$l100\t$l50\t$orconum\t$orconumcom\t$ninec\t$ninei\t$nineic\t$ninep\t$compard\t$ftot\t$fcom\t$fpar\t$fparc\t$fcomp\t$fpse\t$fninec\t$fninei\n";
#print Resultssum "$totalor\t$com\t$comd\t$par\t$pard\t$parc\t$pse\t$psed\t$psec\t$l200\t$l100\t$l50\t$orcount\t$orconum\t$orconumcom\t$ninec\t$ninei\t$nineic\t$ninep\t$sugar\t$bitter\t$co2\t$gr43\t$compard\t$ftot\t$fcomp\t$denovoc\t$fpar\t$denovop\t$fpse\t$fninec\t$fninei\t$frpartial\t$frpse\n";
print Resultssum "$totalor\t$orcount\t$l50\t$ftot\t$fcomp\t$denovoc\t$fpar\t$fpse\t$frpartpse\n";
close Resultssum;

