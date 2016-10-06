
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

my @bdrs = split /-/, $region;

print ("Beginning processing> ".$region." ".$file."\n");

my $filename = $Rof;
open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";

my $fileseq = $Sof;
open(my $fh2, '>', $fileseq) or die "Could not open file '$filename' $!";

my @frequencies;

open my $in, "<:encoding(utf8)", $file or die "$file: $!";
while (my $line = <$in>) {
    
    my @tempsplit = split /\t/, $line;

    if (length($tempsplit[9]) > $bdrs[0]){
	my $sequence = substr($tempsplit[9], $bdrs[0], ($bdrs[1]-$bdrs[0]));	
	print ($sequence."\n");
	my $sout = ">".$ident."\n".$sequence."\n";
	print $fh2 $sout;
    }

    if ($ident == $tempsplit[2]){
	
	if ($tempsplit[5] !~ /^S/){


	    (my @Ms) = $tempsplit[5] =~ /\d{1,3}M/g;
	    (my @Is) = $tempsplit[5] =~ /\d{1,3}[I]/g;
	    (my @Ds) = $tempsplit[5] =~ /\d{1,3}[D]/g;
	    (my @all) = $tempsplit[5] =~ /\d{1,3}/g;
	    
### za prvi entry pogleda, kje je glede na region. ce je pod, potem bo prvemu pristet drugi, dokler ni nad prvo region. Pogledamo koliko nad region je, pristejemo to k dolzini	    
	    
	    my $Msum = 0;
	    my $baseLen = 0;
	    my $insertion = 0;
	    my $deletion = 0;
	    foreach $entry (@all){
		foreach $match (@Ms){
		    if ($entry == $match){
			print($entry." ".$match."\n");
			$Msum += $match;
			#print($Msum."\n");
			if ($bdrs[0] < $Msum && $bdrs[1] > $Msum){
			    print ("lower bound: ".$bdrs[0]." ".$Msum."\n");
			    $baseLen += ($Msum - $bdrs[0]);
			}elsif ($Msum < $bdrs[1] && $Msum > $bdrs[0]){
			    print ("MIddle bound: ".$bdrs[1]." ".$Msum."\n");
			    $baseLen += ($match-$bdrs[0]);
			}elsif ($Msum > $bdrs[1]){
			    print ("Higher bound..\n");
			    $baseLen = ($bdrs[1] - ($bdrs[0]));
			}		   
		    }
		    
		}	   
	    }
	    if ($Msum >= $bdrs[1]){
		print ("....BEG....\n");
		#print ($tempsplit[5]."\n");
		for $sig (@Is){
		    $sig =~ s/[A-Z]+//;
		    print ($sig."iiiiii\n");
		    $baseLen += $sig;
		}		# 
		for $dig (@Ds){
		    $dig =~ s/[A-Z]+//;
		    print ($dig."ddddd\n");
		    $baseLen -= $dig;
		}

		push(@frequencies, $baseLen.","); 
		print ("\n base length: ".$baseLen."\n");
		print ("Match sum: ".$Msum."\n");	
		print ("....END....\n");    
	    }	    
	}
    }
}

close $in;

#close $fhR;
#close $fhS;
#print (@frequencies);
# $lowerBound = $bdrs[1]-$bdrs[0];
# foreach my $x ($lowerBound..($lowerBound+20)){
#     print $x.",";  
#     my $count = grep ($x, @frequencies);
#     print $count.",";
 
# }
    
# my $count = ($string =~ s/a/a/g);

# foreach my $it (@frequencies){
#     print $fh $it;
# }


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
    print $fh "$key:$value\n";
}
close $fh;
close $fh2;
print "done writing to files.\n";
