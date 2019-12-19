#! /usr/bin/env perl

##script to create YAML for Varlociraptor based on contents of dir
use strict;
use warnings;
use Cwd;

##inputs is:
my $GERMID=$ARGV[0];

##define output, concatenate per sample.observations.vcf encountered starting
my $OUTPUT="samples:\n";
my $OUTCALL="varlociraptor call variants --output cons.calls.bcf generic --scenario varlociraptor_call.yaml --obs $GERMID=$GERMID" . ".observations.vcf ";

##standard entry for tumour samples
##NB YAML takes double space instead of tabs
my $STDT="    resolution: 1\n    universe: \"[0.0,1.0]\"\n    contamination:\n      by: $GERMID\n      fraction: 0.25\n";

##put in germline
$OUTPUT.="  $GERMID:\n    resolution: 1\n    universe: \"0.0 | 0.5 | 1.0 | ]0.0,0.5[\"\n";

##iterate over vcfs, catch basename (sampleID) and write entry
my @HOLDIDs;
my @VCFs=glob("*.vcf");

foreach my $vcf (@VCFs){
  my @s=split(/\./,$vcf);
  if($s[0] ne $GERMID){
    $OUTPUT.="  " . $s[0] . ":\n$STDT";
    push(@HOLDIDs, $s[0]);
    $OUTCALL.=$s[0] . "=" . $vcf . " ";
  }
}

##set entries to test; first germline, somatic normal
$OUTPUT.="\nevents:\n  germline:  \"$GERMID:0.5 | $GERMID:1.0\"\n  somatic_normal:  \"$GERMID:]0.0,0.5[\"\n";

foreach my $id (@HOLDIDs){
  $OUTPUT.="  somatic_" . $id . ":  \"$GERMID:0.0 & " . $id . ":]0.0,1.0]\"\n";
}

$OUTCALL.="\n";
open(OUTCAL,">varlociraptor_scenario_call.sh");
print OUTCAL $OUTCALL;
close OUTCAL;
print $OUTPUT;
