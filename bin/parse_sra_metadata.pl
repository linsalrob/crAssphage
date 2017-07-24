#!/usr/bin/perl -w
#

use strict;

open(IN, "biosample_set.xml") || die "can't open biosample_set.xml";
my $id=""; my $sra=""; my $country="";
while (<IN>) {
	next if (index($_, "BioSampleSet") >= 0);
	if (index($_, "</BioSample>") >= 0) {
		if ($country) {
			print join("\t", $id, $sra, $country), "\n";
		}
		$id=""; $country=""; $sra="";
	}

	if (index($_, "<BioSample") >= 0) {
		# <BioSample submission_date="2008-04-04T08:44:25.263" last_update="2015-03-20T11:30:45.343" publication_date="2008-04-04T08:44:25.357" access="public" id="3" accession="SAMN00000003">
		m/accession="(.*?)"/;
		$id = $1;
		if (!$id) {
			print STDERR "ACCESSION: $_";
		}
	}


	if (index($_, '<Id db="SRA"') >= 0) {
		# <Id db="SRA">SRS000002</Id>
		if (m#>(.*?)</#) {
			$sra=$1;
		}
		else {
			print STDERR "SRA: $_";
		}
	}

	if (m/country/i) {
		if (m#<Attribute#) {
			chomp;
			s/^\s*<Attribute.*?>//;
			s#</Attribute>##;
			$country=$_;
		}
		elsif (m#<Paragraph#) {
			if (m/country\s*:([:\w\s]+)/i) {$country = $1}
			elsif (m/country\s*=([:\w\s]+)/i) {$country = $1}
			else {
				print STDERR "PARAGRAPH: $_";
			}
		}
		else {
			print STDERR "OTHER: $_";
		}
	}
}


