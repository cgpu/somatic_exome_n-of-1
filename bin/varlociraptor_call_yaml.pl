#! /usr/bin/env perl

##script to create YAML for Varlociraptor based on contents of dir
use strict;
use warnings;
use Cwd;

##inputs is:
my $GERMID=$ARGV[0];

##define output, concatenate per sample.observations.vcf encountered starting
my $OUTPUT="samples:\n";
my $OUTCALL="varlociraptor call variants generic --scenario varlociraptor_call.yaml --obs germ=" . $GERMID . ".observations.vcf";

##standard entry for tumour samples
my $STDT="\t\tresolution: 10\n\t\tuniverse: \"[0.0,1.0]\"\n\t\tcontamination:\n\t\t\tby: normal\n\t\t\tfraction: 0.25\n";

##put in germline
$OUTPUT.="\tgerm:\n\t\tresolution: 10\n\t\tuniverse: \"0.0 | 0.5 | 1.0 | ]0.0,0.5[\"\n";

##iterate over vcfs, catch basename (sampleID) and write entry
my @HOLDIDs;
my @VCFs=glob("*.vcf");

foreach my $vcf (@VCFs){
  my @s=split(/\./,$vcf);
  $OUTPUT.="\t" . $s . ":\n$STDT";
  push(@HOLDIDs, $s);
  $OUTCALL.=$s . "=" . $vcf . " ";
}

##set entries to test; first germline, somatic normal
$OUTPUT.="\nevents:\n\tgermline:\t\"germ:0.5 | germ:1.0\n\tsomatic_normal\tgerm:\"]0.0,0.5[\"\n";

foreach my $id (@HOLDIDs){
  $OUTPUT.="\tsomatic_" . $id . ":\t\"germ:0.0 & " . $id . ":]0.0,1.0]\"\n";
}

$OUTCALL.=" > cons.calls.vcf";
open(OUTCAL,">varlociraptor_scenario_call.sh");
print OUTCAL $OUTCALL;
close OUTCAL;
print $OUTPUT;
