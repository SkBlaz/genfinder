
## this perl script enables instant counting of pattern sequences amongst
## BAM files..


#najprej converta v SAM za boljse procesiranje 

#samtools view -h -o out.sam in.bam
#mawk '{print $3,"  ",$6, "   ",$10}' ms2.sam 
# hashmap z encodanimi opcijami kako in kaj narediti!

#$ARGV[0]



# input parametri za ta program bodo npr. velikost regije, ime markerja, output fajl.
# input je identifier regija ngs.sam

my $ident = $ARGV[0];
my $region = $ARGV[1];
my $file = $ARGV[2];
my $Rof = $ARGV[3]; #region outfile
my $Sof = $ARGV[4]; #sequence outfile
my $refLen = $ARGV[5]; #sequence outfile

my @bdrs = split /-/, $region;

print ("Beginning processing> ".$region." ".$file."\n");

my $filename = $Rof;
open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";

my $fileseq = $Sof;
open(my $fh2, '>', $fileseq) or die "Could not open file '$filename' $!";

my @frequencies;

open my $in, "<:encoding(utf8)", $file or die "$file: $!";

my $resultcount = 0;
while (my $line = <$in>) {
    
    my @tempsplit = split /\t/, $line;

    ## ce je ident enak imenu npr VVS2
    #dodati iskanje referencnih sekvenc.
    #my $refLen = 137;
    if ($ident eq $tempsplit[2]){
	my $threshold = 0;
	my $orientation = "";
	#print (length($tempsplit[9])." ".$tempsplit[1]."\n");
	if ($tempsplit[1] == 16){
	    if (($refLen-length($tempsplit[9])) <= $bdrs[0]){
		$threshold = 1;
		$orientation = "rev";
		#$resultcount += 1;
	    }
	}else{
	    if (length($tempsplit[9]) >= $bdrs[1]){	
		$threshold = 1;
		$orientation = "for";
		#$resultcount += 1;
	    }
	}
	if ($threshold==1){	
	    #if ($tempsplit[5] !~ /(.|..)S/){

# 		if ($orientation == "for"){
# 		    (my @Ms) = $tempsplit[5] =~ /\d{1,3}M/g;
# 		    (my @Is) = $tempsplit[5] =~ /\d{1,3}[I]/g;
# 		    (my @Ds) = $tempsplit[5] =~ /\d{1,3}[D]/g;
# 		    (my @all) = $tempsplit[5] =~ /\d{1,3}/g;
		    
# 		    ### za prvi entry pogleda, kje je glede na region. ce je pod, potem bo prvemu pristet drugi, dokler ni nad prvo region. Pogledamo koliko nad region je, pristejemo to k dolzini	    
		    
# 		    my $Msum = 0;
# 		    my $baseLen = 0;
# 		    my $insertion = 0;
# 		    my $deletion = 0;
# 		    foreach $entry (@all){
# 			foreach $match (@Ms){
# 			    if ($entry == $match){
# 				#print($entry." ".$match."\n");
# 				$Msum += $match;
# 				#print($Msum."\n");
# 				if ($bdrs[0] < $Msum && $bdrs[1] > $Msum){
# 				    #print ("lower bound: ".$bdrs[0]." ".$Msum."\n");
# 				    $baseLen += ($Msum - $bdrs[0]);
# 				}elsif ($Msum < $bdrs[1] && $Msum > $bdrs[0]){
# 				    #print ("MIddle bound: ".$bdrs[1]." ".$Msum."\n");
# 				    $baseLen += ($match-$bdrs[0]);
# 				}elsif ($Msum > $bdrs[1]){
# 				    #print ("Higher bound..\n");
# 				    $baseLen = ($bdrs[1] - ($bdrs[0]));
# 				}		   
# 			    }
			    
# 			}	   
# 		    }
# 		    if ($Msum >= $bdrs[1]){
# 			#print ("....BEG....\n");
# 			#print ($tempsplit[5]."\n");
# 			for $sig (@Is){
# 			    $sig =~ s/[A-Z]+//;
# 			    #print ($sig."iiiiii\n");
# 			    $baseLen += $sig;
# 			}		# 
# 			for $dig (@Ds){
# 			    $dig =~ s/[A-Z]+//;
# 			    #print ($dig."ddddd\n");
# 			    $baseLen -= $dig;
# 			}

# 			my $sequence = substr($tempsplit[9], $bdrs[0], ($bdrs[1]-$bdrs[0]));	
# 			my $sout = ">".$ident."\n".$sequence."\n";
# 			print $fh2 $sout;       	    
# 			push(@frequencies, $baseLen.","); 
# #			$resultcount += 1;
			
# 		    }
# 		}

		if ($orientation == "rev"){
		    (my @Ms) = $tempsplit[5] =~ /\d{1,3}M/g;
		    (my @Is) = $tempsplit[5] =~ /\d{1,3}[I]/g;
		    (my @Ds) = $tempsplit[5] =~ /\d{1,3}[D]/g;
		    (my @all) = $tempsplit[5] =~ /\d{1,3}/g;
		    
		    ### za prvi entry pogleda, kje je glede na region. ce je pod, potem bo prvemu pristet drugi, dokler ni nad prvo region. Pogledamo koliko nad region je, pristejemo to k dolzini	    
		    my $startlen = $refLen-length($tempsplit[9]);
		    my $Msum = $startlen;
		    my $baseLen = 0;
		    my $insertion = 0;
		    my $deletion = 0;
		    foreach $entry (@all){
			foreach $match (@Ms){
			    if ($entry == $match){
				#print($entry." ".$match."\n");
				$Msum += $match;
				#print($Msum."\n");
				if ($bdrs[0] < $Msum && $bdrs[1] > $Msum){
				    #print ("lower bound: ".$bdrs[0]." ".$Msum."\n");
				    $baseLen += ($Msum - $bdrs[0]);
				}elsif ($Msum < $bdrs[1] && $Msum > $bdrs[0]){
				    #print ("MIddle bound: ".$bdrs[1]." ".$Msum."\n");
				    $baseLen += ($match-$bdrs[0]);
				}elsif ($Msum > $bdrs[1]){
				    #print ("Higher bound..\n");
				    $baseLen = ($bdrs[1] - ($bdrs[0]));
				}		   
			    }
			    
			}	   
		    }
		    if ($Msum >= $bdrs[1]){
			#print ("....BEG....\n");
			#print ($tempsplit[5]."\n");
			for $sig (@Is){
			    $sig =~ s/[A-Z]+//;
			    #print ($sig."iiiiii\n");
			    $baseLen += $sig;
			}		# 
			for $dig (@Ds){
			    $dig =~ s/[A-Z]+//;
			    #print ($dig."ddddd\n");
			    $baseLen -= $dig;
			}

			my $sequence = substr($tempsplit[9], $bdrs[0], ($bdrs[1]-$bdrs[0]));	
			my $sout = ">".$ident."\n".$sequence."\n";
			print $fh2 $sout;       	    
			$resultcount += 1;
			push(@frequencies, $baseLen.","); 
		#	$resultcount += 1;
			
		    }		    #completely different processing scheme.

		}
		
	   # }
	    
	}
	
    }
}
close $in;

my %count;
$count{$_}++ foreach @frequencies;

#removing the lonely strings
while (my ($key, $value) = each(%count)) {
    if ($value == 0) {
        delete($count{$key});
    }
}

#output the counts
while (my ($key, $value) = each(%count)) {
    if ($key > 0){
	print $fh  "$key$value\n"; #$fh
    }
}

close $fh;
close $fh2;
print "done writing to files.\n".$resultcount." results found.";

