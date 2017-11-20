#!/usr/bin/perl
use strict;
use warnings;

#Kommandozeilen-Kommunikation
my $console_uniprot;			#Pfad zur UniProt xml Datei

if(@ARGV!=1){
	print "Error! Keine, zuviele oder falsche Parameter Eingabe!" . "\n";
	print "Parameter: " . "\n";
	print "<Path2proteome.swissprot.xml>" . "\n";
	exit;
}
else{
	$console_uniprot = $ARGV[0];
}

#Variablen-Initialisierung
my $row;
my $counter;
my $last_id;
my $last_name;
my $last_species;
my @locations_swissprot = ();
my %evidences;

#-------------------------------------------------------#
#Ausführung der einzelnen Subroutinen (Main-Runner Code)
#-------------------------------------------------------#
print "uniprot_id" . "\t" . "species" . "\t" .  "gene_name" . "\t" . "annotation_source" . "\t" . "subcellular_location" . "\t" . "annotation_evidence" . "\n";

#Einlesen der UniProt xml Datei
my $switcher = 0;
my $switcher_uni = 0;
my $switcher_go = 0;
open(my $input, "<", $console_uniprot) || die "UniProt Entries xml Datei nicht gefunden: $!";
while($row = <$input>){
	if($row =~ "\n"){
		chop $row;
	}

	if($row =~ "<accession>" && $switcher == 0){
		$row =~ /.*<accession>([A-Za-z0-9]+)<\/accession>/;
		if(defined($1)){
			$last_name = ".";
            $last_species = ".";
			$last_id = $1;
			$counter++;
            @locations_swissprot = ();
            %evidences = ();
		}
		else{
			print "ERROR: ID nicht gefunden." . "\n";
			print "Letzte eingelesene Zeile:" . "\n";
			print $row . "\n";
			exit;
		}
		$switcher = 1;
	}
	elsif($row =~ /<name type="primary".*>(.*)<\/name>/){
		$last_name = $1;
	}
    elsif($row =~ /<name type="scientific".*>(.*)<\/name>/){
        $last_species = $1;
    }
	elsif($row =~ /.*<subcellularLocation>/){
		$switcher_uni = 1;
	}
    # ignore locations without evidence
	elsif($row =~ /.*<location evidence="(\d+)".*>(.*)<\/location>/ && $switcher_uni == 1){
        push(@locations_swissprot, [$2, $1]);
		$switcher_uni = 0;
	}
    ## ignore GO terms atm. 
	# elsif($row =~ /.*<dbReference type="GO" id=".*">/){
	# 	$switcher_go = 1;
	# }
	# elsif($row =~ /.*<property type="term" value="C:(.*)"\/>/ && $switcher_go == 1){
	# 	print $last_id . "\t" . $last_species . "\t" . $last_name . "\t" . "GO" . "\t" . $1 . "\n";
	# 	$switcher_go = 0;
	# }
    elsif($row =~ /<evidence key="(\d+)" type="(.+)".*>/) {
        # print $1 . "\t" .  $2;
        $evidences{$1} = $2;
    }
	elsif($row =~ /\<\/entry.*/ && $switcher == 1) {
        foreach my $loc (@locations_swissprot) {
            my $tmp_loc = $loc->[0];
            my $tmp_evidence = $evidences{$loc->[1]};
            print $last_id . "\t" . $last_species . "\t" . $last_name . "\t" . "SwissProt" . "\t" . $tmp_loc . "\t" . $tmp_evidence . "\n";
        }
		$switcher = 0;
	}
}
close($input);


