#usr/bin/perl -w

#********************************************************************************************************************************************************************************************************
# in the v2 version HLA bounderies have been corrected to a subregion
# v3 removed HLA filter entirely.  It still filters for PAR regions and with QuantiSNP for nb of probes, length and score.
#
#*******************This prog was adapted by Maude Auger for internal need.  This file is independent from the rest of CNVision package, it runs by itself and won't call PennCNV, QuantiSNP or java files
#
# --PNformat /media/maude/MaudeLabBackup/labo/myPipeline/CCP_2019v1/testCNVision/RawPennCNVoutput/PC_allCNV.txt /media/maude/MaudeLabBackup/labo/myPipeline/CCP_2019v1/testCNVision/test01 Penn-tag
# --QTformat /media/maude/MaudeLabBackup/labo/myPipeline/CCP_2019v1/testCNVision/QS_allCNV.txt /media/maude/MaudeLabBackup/labo/myPipeline/CCP_2019v1/testCNVision/test01 QS-tag
# --merge /media/maude/MaudeLabBackup/labo/myPipeline/CCP_2019v1/testCNVision/test01/PC_testclean_Formated.txt /media/maude/MaudeLabBackup/labo/myPipeline/CCP_2019v1/testCNVision/test01/QS_testclean_Formated.txt /media/maude/MaudeLabBackup/labo/myPipeline/CCP_2019v1/testCNVision/test01/ merge-tag 
#
#
#  /home/maudonna/projects/rrg-jacquese/All_user_common_folder/RAW_DATA/AutismGenomeProject-AGP/AGP_analysis/v3/CNVision
#  /RQexec/guillaf/GRIP/grip2analysis
#**********************************************************************************************************************************************************************************************************


use Getopt::Long;

#Checks the time
$timeanddate = &showtime;
#This prevents Excel interpreting the date:
#$timeanddate = '\''.$timeanddate.'\'';
#print "Time: $timeanddate\n";


GetOptions(  'pnformat' => \$PNformat,
             'qtformat' => \$QTformat,
             'merge' => \$merge_prog, 
             'help' => \$help,            
);
	

print "=============================================================================\n";
print "   CNV Pipeline by Stephan Sanders        email: stephan.sanders\@yale.edu\n";
print "   and Christopher Mason                  \n";
print "          \n";
print "   reduced to pennCNV formating, QuantiSNP formating and merge by Maude Auger\n";
print "=============================================================================\n";


#**********************************************************************************************************************************************************************************************************************************************************
#**********************************************************************************************************************************************************************************************************************************************************
#**********************************************************************************************************************************************************************************************************************************************************
#**********************************************************************************************************************************************************************************************************************************************************
#                                 
#  PENN CNV FORMAT     ---    PENN CNV FORMAT   ---    PENN CNV FORMAT     ---    PENN CNV FORMAT   ---   PENN CNV FORMAT     ---    PENN CNV FORMAT   ---   PENN CNV FORMAT     ---    PENN CNV FORMAT   ---
# 
#Formats PennCNV files into CNVision input format - original script was not outputing what was expected, therefore the section was completely modified. we added filter to remove HLA and PAR regions.
#
#**********************************************************************************************************************************************************************************************************************************************************
#**********************************************************************************************************************************************************************************************************************************************************
#**********************************************************************************************************************************************************************************************************************************************************
#*********************************************************************************************************************************************************************************************************************************************************

if ($PNformat == 1){
	print "CNVision: Reformatting PennCNV results\n";
	
 
	if (scalar @ARGV != 3) {
		print "PNformat option uses 3 arguments: the path of PennCNV file, the output directory and a cohort identifier that will be in the name of the output file\n";
	} 
	
	$tag = $ARGV[2];
	$pro_dir = $ARGV[1];
	if (! -d $pro_dir) {
   	print "The second argument should be the path of the output directory, you did not give the path for a directory\n";
   	exit;
	}
	$pro_dir = $pro_dir . "/";
	$pnres = $ARGV[0];
	
	
	$PNform = "$pro_dir"."PC_" . $tag . "_CNVisionFormated.txt";
	
	
	open (PNRES, ">$PNform" || die "Could not open $PNform; may already be open\n");
	print PNRES "FID\tSampleID\tChr\tStart\tStop\tType\tConf\tNbProbes\tSize\tStartSNP\tStopSNP\n";
	
	open (PRES, "<$pnres" || die "ERROR: No results file for PennCNV: looking for $pnres\n");
	
	
	while (<PRES>){
		push (@PN, $_);
	}
	close PRES;
	
	
	$theM = 'M';
	$theX = 'X';
	$thechr = 'chr';
	$callcount = 0;
	$xcallcount = 0;
	$excluded = 0;
	@PN = sort { (split /\s+/, $a)[4] cmp (split /\s+/, $b)[4] } @PN;
	
	foreach $line (@PN){
		$linecount++;
		@tempa = split(/\s+/,$line);
		@nextline = split(/\s+/,@PN[$linecount]);
		chomp($tempa[0]);
		@chromoall = split(":", $tempa[0]);
		$chromo = $chromoall[0];
		#print "$chromo\n";
		if ($chromo =~ /$theX/){
			$xcallcount++;
		}
		else{
			$callcount++;
		}
		$chromo =~ s/MT/$theM/gi;
		$chromo =~ s/$thechr//gi;
		$chromo = "chr".$chromo;
		($CNVstart, $CNVstop) = split("-",  $chromoall[1]);
		$CNVstart =~ s/,//gi;
		$CNVstop =~ s/,//gi;
		$tempa[1] =~ s/numsnp=//i;
		$tempa[2] =~ s/length=//i;
		$tempa[2] =~ s/,//gi;
		$tempa[3] =~ s/state5,cn=//i;
		$tempa[3] =~ s/state6,cn=//i;
		$tempa[3] =~ s/state2,cn=//i;
		$tempa[3] =~ s/state1,cn=//i;
		$tempa[3] =~ s/state3,cn=//i;
		$thisone = $tempa[4];
		
#*&*&*&*&*&*&*&*&*&*&*&*&*&* from the FR file name, only keep the ID of the sample *&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&		
#*&*&*&*&*&*&*&*&*&*&*&*&*&* adjust for file extension   *&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&	
		#@tab = split (/_/,$tempa[4]);
		
		$samp_remove1 = ".Signal.txt";
		#$samp_remove1 = ".COPYNUM";
		#print "$samp_remove1\n";
		
		
		
		$tempa[4] =~ s/$samp_remove1//gi;
		#print "$tempa[4]\n";
		
#*&*&*&*&*&*&*&*&*&*&*&*&* adjust for prefix if there is one in the file name: ex, here the prefix is hg19_   *&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&	
		#$samp_removePrefix = "hg19_";		
		#$tempa[4] =~ s/$samp_removePrefix//gi;		


#*&*&*&*&*&*&*&*&*&*&*&*&* adjust for AGP: no PN_ in the FR path *&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&		
		#@tab = split (/PN_/,$tempa[4]);
		@tab = split (/\//,$tempa[4]);
		$tempa[4] = $tab[-1];
		#print "$tempa[4]\n";		
		
		$sample{$tempa[4]} = 0;
		$tempa[5] =~ s/startsnp=//i;
		$tempa[6] =~ s/endsnp=//i;
		$tempa[7] =~ s/conf=//i;
		
		
		#coordonates determined by remapping from hg18 (from Guillaume Huguet original PAR script)
		#PAR1 in hg19 chrX:0-2,699,520
		#PAR2 in hg19 chrX:154,931,042-end
		#PAR3 in hg19 chrX:88,456,803-92,375,509
		#exclude CNV if overlapping a PAR region
		#excluded, 1 is True, 0 is False
		if ($chromo eq "chrX" and $tempa[3] == 2){
			$excluded = &checkPARregion($CNVstart, $CNVstop);
		}
		
		#if not excluded save line
		if ($excluded == 1){
			$excluded = 0;
		}else{
			$tempb = "$tempa[4]\t$tempa[4]\t$chromo\t$CNVstart\t$CNVstop\t$tempa[3]\t$tempa[7]\t$tempa[1]\t$tempa[2]\t$tempa[5]\t$tempa[6]\n";
			push (@PNF, $tempb);
		}
		
		if ($thisone ne $nextline[4]){
			$pncalls{$tempa[4]} = "$callcount\t$xcallcount";
			#print "Sample: $thisone, Next: $nextline[4], Sample:$tempa[4], Calls:$callcount\n";
			$callcount = 0;
			$xcallcount = 0;
		}
		
	}
	
	print PNRES @PNF;
	
	@samplecount = keys %sample;
	$samplecount = @samplecount;
	if ($samplecount == 0){
		print "ERROR: No samples found in PennCNV results file\n";
	}
	
	print "PennCNV reformatted for $samplecount samples\n";
}



#**********************************************************************************************************************************************************************************************************************************************************
#**********************************************************************************************************************************************************************************************************************************************************
#**********************************************************************************************************************************************************************************************************************************************************
#**********************************************************************************************************************************************************************************************************************************************************
#                                 
#  QUANTI SNP FORMAT     ---    QUANTI SNP FORMAT   ---    QUANTI SNP FORMAT     ---    QUANTI SNP FORMAT   ---   QUANTI SNP FORMAT     ---    QUANTI SNP FORMAT   ---   QUANTI SNP FORMAT     ---    QUANTI SNP FORMAT   ---
#
#Formats QuantiSNP files into CNVision input format - this is a completely remodeled QTformat section. we added filter to remove HLA and PAR regions and quality filters.
#
#**********************************************************************************************************************************************************************************************************************************************************
#**********************************************************************************************************************************************************************************************************************************************************
#*********************************************************************************************************************************************************************************************************************************************************
#Formats QuantiSNP v1.1 files into Merge format
if ($QTformat == 1){
	print "CNVision: Reformatting Quanti SNP results\n";
	
	if (@ARGV != 3) {
		print "QTformat option uses 3 arguments: the path of QuantiSNP file, the output directory and a cohort identifier that will be in the name of the output file\n";
	} 
	
	$tag = $ARGV[2];
	$pro_dir = $ARGV[1];
	if (! -d $pro_dir) {
   	print "The second argument should be the path of the output directory, you did not give the path for a directory\n";
   	exit;
	}
	$pro_dir = $pro_dir . "/";
	$pnres = $ARGV[0];
	
	
	#quality filters
	$minLength = 1000;
	$minProbes = 3;
	$minScore = 15.0;
	
	

	$PNform = "$pro_dir"."QS_" . $tag . "_CNVisionFormated.txt";

	
	open (PNRES, ">$PNform" || die "Could not open $PNform in QSformat; may already be open\n");
	print PNRES "FID\tSampleID\tChr\tStart\tStop\tType\tConf\tNbProbes\tSize\tStartSNP\tStopSNP\n";
	
	open (PRES, "<$pnres" || die "ERROR: No results file for QuantiSNP: looking for $pnres\n");
	
	
	while (<PRES>){
		push (@PN, $_);
	}
	close PRES;
	
	$theM = 'M';
	$theX = 'X';
	$thechr = 'chr';
	$callcount = 0;
	$xcallcount = 0;
	$excluded = 0;
	$header = 1;
	
	
	foreach $line (@PN){
		if ($header == 1){
			$header = 0;
			next;
					
		}
		$linecount++;
		@tempa = split(/\t/,$line);
		@nextline = split(/\t/,@PN[$linecount]);
		
		$ID = $tempa[0];
		$sample{$ID} = 0;
		$chromo = $tempa[1];
		#print "$chromo\n";
		if ($chromo =~ /23/){
			$chromo = $theX;
			$xcallcount++;
		}
		else{
			$callcount++;
		}
		$chromo =~ s/MT/$theM/gi;
		$chromo =~ s/$thechr//gi;
		$chromo = "chr".$chromo;
		$tempa[2] = $tempa[2] + 0;
		$CNVstart = $tempa[2];
		$tempa[3] = $tempa[3] + 0;
		$CNVstop = $tempa[3];
		$StartSNP = $tempa[4];
		$StopSNP = $tempa[5];
		$tempa[6] = $tempa[6] + 0;
		$Size = $tempa[6];
		if ( $Size < $minLength ){
			$excluded = 1;	
		}
		$tempa[7] = $tempa[7] + 0;
		$NbProbes = $tempa[7];
		if ( $NbProbes < $minProbes ){
			$excluded = 1;	
		}
		$tempa[8] = $tempa[8] + 0;
		$Type = $tempa[8];
		$tempa[2] = $tempa[2] + 0.00;
		$Score =$tempa[9];
		if ( $Score < $minScore ){
			$excluded = 1;
		}
		
		
		#if ($linecount < 20){
		#	$sizeCond = $Size < $minLength;
		#	$probeCond = $NbProbes < $minProbes;
		#	$scoreCond = $Score < $minScore;	
		#	print "size: ($Size; $sizeCond), probres: ($NbProbes; $probeCond), score: ($Score, $scoreCond)\n";
		#}
		
		
		if ($excluded == 0){
			#coordonates determined by remapping from hg18 (from Guillaume Huguet original PAR script)
			#PAR1 in hg19 chrX:0-2,699,520
			#PAR2 in hg19 chrX:154,931,042-end
			#PAR3 in hg19 chrX:88,456,803-92,375,509
			#exclude CNV if overlapping a PAR region
			#excluded, 1 is True, 0 is False
			if ($chromo eq "chrX" and $Type == 2){
				$excluded = &checkPARregion($CNVstart, $CNVstop);
			}
		}
		
		
		
		if ($excluded == 1){
			$excluded = 0;
		}else{
			$tempb = "$ID\t$ID\t$chromo\t$CNVstart\t$CNVstop\t$Type\t$Score\t$NbProbes\t$Size\t$StartSNP\t$StopSNP\n";
			push (@PNF, $tempb);
		}
		
		if ($thisone ne $nextline[0]){
			$pncalls{$tempa[0]} = "$callcount\t$xcallcount";
			#print "Sample: $thisone, Next: $nextline[4], Sample:$tempa[4], Calls:$callcount\n";
			$callcount = 0;
			$xcallcount = 0;
		}
		
	}
	
	print PNRES @PNF;

	
	@samplecount = keys %sample;
	$samplecount = @samplecount;
	if ($samplecount == 0){
		print "ERROR: No samples found in Quanti SNP results file\n";

	}
	print "QUantiSNP reformatted for $samplecount samples\n";
}


#**********************************************************************************************************************************************************************************************************************************************
#**********************************************************************************************************************************************************************************************************************************************
# new addition : identifies PAR1, PAR2 and PAR3 region, in order to remove them from the results and store them in a separate file.
#**********************************************************************************************************************************************************************************************************************************************
#***********************************************************************************************************************************************************************************************************************************************
#By remapping from hg18 (from Guillaume Huguet original PAR script)
#PAR1 in hg19 chrX:0-2,699,520
#PAR2 in hg19 chrX:154,931,042-end
#PAR3 in hg19 chrX:88,456,803-92,375,509
#exclude CNV if overlapping a PAR region
#excluded, 1 is True, 0 is False
sub checkPARregion{
	$start = $_[0];
	$stop = $_[1];
	$exclusion = 0;
	
	if ( $start <= 2699520 ){
		# region PAR1
		$exclusion = 1;	
	}elsif ( $stop >= 154931042 ){
		# region PAR2
		$exclusion = 1;
	}else{
		# region PAR3
		if ( $start >= 88456803 and $start <= 92375509 ){
			$exclusion = 1;	
		}elsif ( $stop >= 88456803 and $stop <= 92375509 ){
			$exclusion = 1;
		}elsif ( $start < 88456803 and $stop > 92375509 ){
			$exclusion = 1;
		}
		
		# si on se retrouve ici c'est que le CNV n'est pas dans une région PAR
	}
	return $exclusion;
}



#Merges multiple files of regions 
if ($merge_prog == 1){
	print "CNVision: merging files\n";

	if (@ARGV != 4) {
		print "merge option requires 4 arguments:\n";
		print "1- The path of PennCNV file formated for CNVision (--PNformat)\n";
		print "2- The path of QuantiSNP file formated for CNVision (--QTformat)\n";
		print "3- The path of the output directory\n";
		print "4- A cohort identifier, that will be in the name of the output file\n";
	} 
	
	$tag = $ARGV[3];
	$pro_dir = $ARGV[2];
	if (! -d $pro_dir) {
   	print "The third argument should be the path of the output directory, you did not give the path for a directory\n";
   	exit;
	}
	$pro_dir = $pro_dir . "/";
	
	$merge_dir = $ARGV[2];
	if (! -d $merge_dir) {
   	print "The third argument should be the path of the output directory, you did not give the path for a directory\n";
   	exit;
	}
	$merge_dir = $merge_dir . "/";
	push (@filesToMerge, $ARGV[0]);
	push (@filesToMerge, $ARGV[1]);

	
	$famID =    0;
	$sampleID = 1;
	$chromo =   2;
	$start =    3;
	$stop =     4;
	$cnvtype =  5;
	$score =    6;
	$numsnps =  7;
	$size =     8;
	$rs1 =      9;#le cn algo1 *****
	$rs2 =     10;#le cn algo2 *****
	
	
	#Initialise UniqueID counter
	$uidc = 0;
	$chrword = "chr";
	
	#open input files 
	#Put the CNVs into an array but add an ID made up of the algo and uniqueID counters	
	# les 2 fichiers seront stocké dans CNVall, mais identifié par les 2 colonnes ajoutées
	foreach $file (@filesToMerge){ 
		print "Loading file $file\n";
		open (INPUT, "<$file") || die "Could not open file $file; check it is spelled correctly\n";
		{ local $/=undef;  $input1=<INPUT>; }
		close INPUT;
		($\) = (/(\r\n|\r|\n)/);
		@CNV=split (/\r\n|\r|\n/, $input1); chomp @CNV;
		chomp;
		
		@fileTab = split (/\//,$file);
		$file = $fileTab[-1];
		$firsttwo = substr($file, 0,2);
print "***************  les fichiers sont $firsttwo  *********************************\n";
		if ($firsttwo{$firsttwo} == 1){
			print "ERROR: The first two letter of the files being merged are the same\n";
		} 
		$firsttwo{$firsttwo} = 1;
		$firstnames = $firstnames.$firsttwo."_";
		$outputName = "CNVisionMerged_".$firstnames.$tag.".txt";
		
		#******tout ca ci-bas c'est pour s'assurer de mettre chr devant le chromosome nb, puis ajouter 2 colonnes au début PN00000001 puis PN ou QT00012345 puis QT  *********
		for ($c=0; $c<=$#CNV; $c++){
			if($CNV[$c] =~ /FID/i or $CNV[$c] =~ /IID/i or $CNV[$c] =~ /sample/i){}
			else{
				@tab = split (/\t/,$CNV[$c]);
				$tab[$chromo] =~ s/$chrword//gi;
				$tab[$chromo] = "chr".$tab[$chromo];
				$CNV[$c] = join ("\t",@tab);
				$uidc++; 
				$uidc = sprintf ("%07d", $uidc);# 7 leading zeros max, all nb have 8 digits
				$CNV[$c] = $firsttwo . $uidc . "\t" . $firsttwo . "\t" . $CNV[$c];
				push (@CNVall, $CNV[$c]);
			}
		}
	}
	
	@CNV = ();
	
	#Two columns have been added so numbers are corrected
	$unique   = 0;
	$algo     = 1;
	$famID    = $famID + 2;
	$sampleID = $sampleID + 2;
	$chromo   = $chromo + 2;
	$start    = $start + 2;
	$stop     = $stop + 2;
	$cnvtype  = $cnvtype + 2;
	$score    = $score + 2;
	$numsnps  = $numsnps +2 ;
	$size     = $size + 2;
	
	#I then need to sort the CNVs by sample, chr, start
	print "Sorting CNVs\n";
	
	#Sort by Start (numbers)
	@CNVall = sort { (split '\t', $a)[$start] <=> (split '\t', $b)[$start] } @CNVall;
	
	#Sort by Chromosome (string)
	@CNVall = sort { (split '\t', lc $a)[$chromo] cmp (split '\t', lc $b)[$chromo] } @CNVall;
	
	#Sort by SampleID (string)
	@CNVall = sort { (split '\t', lc $a)[$sampleID] cmp (split '\t', lc $b)[$sampleID] } @CNVall;
		
		
	print "Looking for CNVs that match\n";
	
	$error = 0;
	$track = 0;
	
	#Does each CNV line up with the one below?
	foreach $line (@CNVall){
		@before = split ("\t", $CNVall[$track - 1]);
		@now = split ("\t", $CNVall[$track]);
		@after = split ("\t",$CNVall[$track + 1]);
		$track++;
		$sampcnvcount++;
		
		if (($now[$stop] - $now[$start]) < 0){
			print "ERROR: CNV $now[0] has a stop co-ordinate lower than the start co-ordinate\n";
			$error++;
		}
		#print "\n\nNew CNV:\n";
		#print "Before: @before\n";
		#print "Now: @now\n";
		#print "After: @after\n";
			
		#If this CNV is from a differet sample or chromosome or if it does not overlap with the last CNV or maxbase then start a new set
		if ($now[$sampleID] ne $before[$sampleID]
		or lc $now[$chromo] ne lc $before[$chromo]
		or ($now[$start] > $before[$stop]
		and $now[$start] > $maxbase)){
			$maxbase = 0;
			@match = ();
			#print "New sample, chromosome or set\n";
			push (@match, $line);
		}
		
		#If it is the same sample and chromosome and this CNV overlaps with the last or with maxbase
		elsif ($now[$sampleID] eq $before[$sampleID] 
		and lc $now[$chromo] eq lc $before[$chromo]
		and ($now[$start] <= $before[$stop]
		or $now[$start] <= $maxbase)){
			
			#All matches get put into an array called @match
			push (@match, $line);
			#print "This was a match\n";
		}
		
		#If the CNV does not fit into either of the above groups then there is an error
		else{
			print "Error with the following CNV, its status could not be determined:\n$line\n";
		}
		
		#If this CNV does not overlap with the next CNV then @match is processed
		if ($now[$sampleID] ne $after[$sampleID]
		or lc $now[$chromo] ne lc $after[$chromo]
		or ($now[$stop] < $after[$start]
		and $maxbase < $after[$start])
		or $after[$start] == ""){
			#Process @match
			&setexaminem;
			$locicounter++;
			#print "Last match in a set, dealing with matches\n@match\n";
			}
		
		#The numbers of loci with all, 1, 2 or 3 algos are entered into the corresponding sample name
		if ($now[$sampleID] ne $after[$sampleID]){
			$sampcnv{$now[$sampleID]} = $sampcnvcount;
			$sampcnvcount = 0;
		}
		
		#Makes Maxbase the highest it can be
		if ($now[$stop] > $maxbase){
			$maxbase = $now[$stop];
		}
	}
	
	&smallinquire(@smout_list);
	
	$timeanddate = &showtime;
	$totalsamp = keys %sampcnv;
	
	@samplist = keys %sampcnv;
	
	foreach (@samplist){
		@tab = split (/\t/,$_);
		$allcnv = $sampcnv{$tab[0]};
		$allloci = $totalloci{$tab[0]};
		$xallloci = $xtotalloci{$tab[0]};
		$twoloci = $samploci2{$tab[0]};
		$oneloci = $samploci1{$tab[0]};
		
		
		
		if ($allloci ne ""){
			print SAMP "$_\t$allcnv\t$allloci\t$xallloci\t$threeloci\t$twoloci\t$oneloci\n";
			$samplecount++;
		}
		else {
			print SAMP "$_\t0\t0\t0\t0\n";
			print "ERROR: No CNVs found by any algorithms for $tab[0]\n";
		}
	}
	if ($samplecount == $totalsamp){
		print "Merge complete for $totalsamp samples\n";
	}
	else {
		print "ERROR: Merge complete, however CNVs were present for samples that were not listed in the Samplelist file\n";
	}
	
	sub setexaminem{
		#chomp removes new line characters
		chomp @match;
		#this will give the length of match
		$matches = @match;
		
		#If there is only one CNV it is simply printed out
		if ($matches == 1){
			$now[$chromo] =~ s/chr//i;
			$mout_line = "$now[$famID]\t$now[$sampleID]\tchr$now[$chromo]\t$now[$start]\t$now[$stop]\t$now[$cnvtype]\t$now[$score]\t$now[$numsnps]\t$now[$size]\t\t$now[$start]\t$now[$stop]\t$now[$size]\t$now[$cnvtype]\t1\t$now[$algo]\t\t$now[$algo]\t$now[$cnvtype]\t$now[$score]\t$now[$numsnps]\t$now[$size]\t$now[$unique]";
			push(@mout_list, $mout_line);
			@match = ();
		}
		
		%coord = ();
		%locitype = ();
		%sizes = ();
		%snps = ();
		%scores = ();
		$counter = 0;
		
		#Put the start and stop positions into a new array and determine the characteristics of this group of matches
		foreach (@match){
			#$_ is the default variable if values are not assign to a var by default it will be to $_ . every function applied to nothing is applied to $_. $_ doesn't have to be explicitly written.
			#all var of ref: start, stop, etc are indexes of their position in each line that was stored in match
			@temp = split ("\t", $_);
			$temp[$start] =~ s/^\s+//;
			$temp[$stop] =~ s/^\s+//;
			chomp $temp[$start];
			chomp $temp[$stop];
			$coord{$temp[$start]} = 0;
			$coord{$temp[$stop]} = 0;
			$locitype{$temp[$cnvtype]} = 0;
			$sizes{$temp[$size]} = 0;
			$scores{$temp[$score]} = 0;
			$snps{$temp[$numsnps]} = 0;
			#les ID et le chr sont les memes pour tous les matchs, donc on utilise une variable
			$family = $temp[$famID];
			$sample = $temp[$sampleID];
			$chrm = $temp[$chromo];
		}
		
		#remove the prefix: chr
		$chrm =~ s/chr//i;
		
		#Sort and determine the loci start and end
		@coord = keys %coord;
		@coord = sort {$a <=> $b} @coord;
		$locistart = $coord[0];
		$locistop = $coord[$#coord];
		$locisize = $locistop - $locistart;
		
		#Sort and determine the loci types present
		@locitypes = keys %locitype;
		@locitypes = sort {$a <=> $b} @locitypes;
		$locitypes = $types = join (", ", @locitypes);
		
		
		#Sort and determine the largest sized call
		@sizes = keys %sizes;
		@sizes = sort {$a <=> $b} @sizes;
		$maxsize = $sizes[-1];
		
		#Sort and determine the highest score in a call
		@scores = keys %scores;
		@scores = sort {$a <=> $b} @scores;
		$maxscore = $scores[-1];
		
		#Sort and determine the call with the most snps
		@snps = keys %snps;
		@snps = sort {$a <=> $b} @snps;
		$maxsnp = $snps[-1];
		
		#print "Within this set: Start = $locistart, End = $locistop, Types = $locitypes, Maxsize = $maxsize, MaxSNP = $maxsnp\n";
		
		$algocount = 0;
		%algos = ();
		%uid = ();
		%type = ();
		
		############################################################################################################################
		
		#Look for matches at each line of the co-ord
		foreach $bp (@coord){
			$lastone = $coord[$counter - 1];
			$thisone = $coord[$counter];
			$nextone = $coord[$counter + 1];
			
			$counter++;
			
			##########################################################################################################################
			
			#For each coord look for matches in the start and stop coords
			foreach $line4 (@match){
				@cnvm = split ("\t", $line4);
				#If the coord matches the start of a cnv the uid and features are added to hashes as the key and value
				if($thisone == $cnvm[$start]){
					$algocount++;
					chomp $cnvm[$unique];
					$uid{$cnvm[$unique]} = "$cnvm[$algo]\t$cnvm[$cnvtype]\t$cnvm[$score]\t$cnvm[$numsnps]\t$cnvm[$size]\t$cnvm[$unique]";
				}
				
				#If the coord matches the end of a CNV the uid and algo are removed from the hashes
				if ($thisone == $cnvm[$stop]){
					$algocount = $algocount - 1;
					chomp $cnvm[$unique];
					delete ($uid{$cnvm[$unique]});
				}
			}
			
			%type = ();
			%algoname = ();
			@features = ();
			
			@uidkeys = keys %uid;
			
			foreach $key (@uidkeys){
				$temp = $uid{$key};
				chomp $temp;
				push (@features, $temp);
				@bits = split ("\t", $temp);
				$type{$bits[1]} = 1;
				$algoname{$bits[0]} = 1;
			}
			
			#Work out the size of this section
			$musize = $nextone - $thisone;
			
			#Work out the types present
			@type = keys %type; 
			chomp @type;
			@type = sort @type;
			$types = join (", ", @type);
			
			#Work out the algorithms present
			@algoname = keys %algoname;
			@algoname = sort @algoname;
			$algonames = join (", ", @algoname);
			$algonamesize = @algoname;
			if ($algonamesize > 1){
				foreach $key (@uidkeys){
					$uidmatch{$key} = 1;
				}
			}
			
			#Printout this section of the match
			$printout1 = "$family\t$sample\tchr$chrm\t$locistart\t$locistop\t$locitypes\t$maxscore\t$maxsnp\t$locisize\t\t$thisone\t$nextone\t$musize\t$types\t$algocount\t$algonames\t\t";
			chomp @features;
			$printout2 = join ("\t", @features);
			unless ($nextone == ""){
				$new_mout_line = $printout1 . $printout2;
				push(@mout_list, $new_mout_line);
			}
		}
	}
	foreach $line (@CNVall){
		chomp $line;
		$line =~ s/\\n//gi;
		@temp = split ("\t", $line);
		$uidline = "$line\t$uidmatch{$temp[0]}";
		push(@UID_list, $uidline)
	}
}


#Subroutine to give the current time
sub showtime {
	@months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	@weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	$year = 1900 + $yearOffset;
	$hour= sprintf("%02d", $hour);
	$minute= sprintf("%02d", $minute);
	$second= sprintf("%02d", $second);
	$smallTime = "h$hour" . "m$minute";
	$theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
	$timeanddate = "$hour:$minute $dayOfMonth/$months[$month]/$year";
	$stamp = "$months[$month]$dayOfMonth,$year";
	return $timeanddate;
}


#Subroutine used for small inquire
sub smallinquire {
	
	$output = $merge_dir . "Sum_". $outputName;
	
	open (OUT, ">$output") || die "Cannot open file $output, may already be open\n";
	print OUT "FID\tSampleID\tCHR\tSTART\tSTOP\tType\tSCORE\tSNP\tSize\t#Algos\tAlgos\t%Three Algs\tTwoAlgs\t%One Alg\n";

	@CNV = @mout_list;
	$extra = 0;
	
	#les index pour printout1
	$famID = 0 + $extra;
	$sampleID = 1 + $extra;
	$chromo = 2 + $extra;
	$lstart = 3 + $extra;
	$lstop = 4 + $extra;
	$cnvtype= 5 + $extra;
	$maxsize = 6 + $extra;
	$maxsnp = 7 + $extra;
	$lsize = 8 + $extra;
	$start = 10 + $extra;
	$stop = 11 + $extra;
	$size = 12 + $extra;
	$type = 13 + $extra;
	$algonum = 14 + $extra;
	$algos = 15 + $extra;
	
	###############################################################################################
	print "Calculating numbers of loci present\n";
	#print "Sorting input file by Sample, Chromosome and Loci Start\n";
	#Sort by Start (numbers)
	@CNV = sort { (split '\t', $a)[$start] <=> (split '\t', $b)[$start] } @CNV;
	
	#Sort by Chromosome (string)
	@CNV = sort { (split '\t', lc $a)[$chromo] cmp (split '\t', lc $b)[$chromo] } @CNV;
	
	#Sort by SampleID (string)
	@CNV = sort { (split '\t', lc $a)[$sampleID] cmp (split '\t', lc $b)[$sampleID] } @CNV;
	
	#Split CNV list into loci then process each sample as it ends
	#print "CNVs sorted, now examining each loci\n";
	$tracker = 0;
	foreach (@CNV){
		@before = split ("\t", $CNV[$tracker - 1]);
		@now = split ("\t", $CNV[$tracker]);
		@after = split ("\t", $CNV[$tracker + 1]);
		$tracker++;
		
		#print "\nNew CNV:\n";
		#print "Before: @before\n";
		#print "Now: @now\n";
		#print "After: @after\n";
		
		push (@loci, $_);
		#print "Now: $now[$lstart], After: $after[$lstart]\n";
		
		#in the Merged......txt file if consecutive CNV are not the same then
		if ($now[$sampleID] ne $after[$sampleID]
		or $now[$chromo] ne $after[$chromo]
		or $now[$lstart] != $after[$lstart]
		or $after[$sampleID] eq ""){
			#print "Now: $now[$sampleID], After: $after[$sampleID]\n";
			#print "This is the next loci:\n";
			&lociprosi;
			$locicount++;
			@loci = ();
		}
		
		#if we changed individual - store info for the current individual
		if ($now[$sampleID] ne $after[$sampleID]){
			$totalloci{$now[$sampleID]} = $totalloci;
			$xtotalloci{$now[$sampleID]} = $xtotalloci;
			$samploci1{$now[$sampleID]} = $samploci1count;
			$samploci2{$now[$sampleID]} = $samploci2count;
			$samploci1count = 0;
			$samploci2count = 0;
			$totalloci = 0;
			$xtotalloci = 0;
		}
	}
	
	sub lociprosi{
		$locicountx++;
		$two = 0;
		$one = 0;
		%maxalgo = ();
		
		foreach $line (@loci){
			@tabs = split ("\t", $line);
			
			if ($tabs[$algonum] == 3){
				print "ERROR: A problem was encountered in lociprosi subroutine, it recognised 3 algo, the maximum is 2.";
				exit;
			}
			if ($tabs[$algonum] == 2){
				$two = $two + $tabs[$size];
			}
			if ($tabs[$algonum] == 1){
				$one = $one + $tabs[$size];
			}
			$maxalgo{$tabs[$algonum]} = $tabs[$algos];
		}
		
		@tab = split ("\t", $loci[0]);
		
		#Percent rare1, de novo and number of algos
		if ($tab[$lsize] > 0){
			$ptwo = sprintf ("%0d", $two / $tab[$lsize] * 100);
			$pone = sprintf ("%0d", $one / $tab[$lsize] * 100);
		}
		
		@algo = keys %maxalgo;
		@algo = sort {$a <=> $b} @algo;
		$maxal = $algo[-1];
		$maxalgs = $maxalgo{$maxal};
		
		#These counters are used to work out the number of loci with different numbers of algorithms
		if ($tab[$chromo] =~ /x/i){
			$xtotalloci++;
		}
		$totalloci++;
		if ($maxal == 1){
			$samploci1count++;
		}
		if ($maxal == 2){
			$samploci2count++;
		}
		
		
		print OUT "$tab[$famID]\t$tab[$sampleID]\t$tab[$chromo]\t$tab[$lstart]\t$tab[$lstop]\t$tab[$cnvtype]\t$tab[$maxsize]\t$tab[$maxsnp]\t$tab[$lsize]\t$maxal\t$maxalgs\t0%\t$ptwo\t$pone%\n";
		
	}
	print "\n";
}



if ($help == 1){
  print "\nTry entering the following commands:\n";
  print "CNV Pipeline:\n";
  print "--merge -------> Merges files in format Family/Sample/Chr/Start/Stop/Type/Score\n";
  print "                 /SNPs/Size; list files for merging after the command\n";
  print "--PNformat --------> Transforms raw PennCNV output into a CNVision format\n";
  print "--QTformat ------> Transforms raw QuantiSNP output into a CNVision format\n";           
}







