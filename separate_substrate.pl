#!/usr/bin/perl
use strict;
use warnings;
use Math::Trig;
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use lib dirname(abs_path $0);

require 'crossover_subroutines.pl';

#----------------------------------------------------------
# User Inputs - Geometry File, Cluster length

my @first_collection = @{read_file("$ARGV[0]/geometry.in")};
my $length = $ARGV[1];

#----------------------------------------------------------
# Isolating the lattice vectors and finding the magnitude
my ($lattice_vectors_ref, $first_collection_ref)  = extract_lattice_vectors(\@first_collection);
my @lattice_vectors = @{$lattice_vectors_ref};
my @magnitudes = @{get_magnitude_of_lattice_vectors(\@lattice_vectors)};

#----------------------------------------------------------
# Generating a rough approximation for the cluster

@first_collection = @{$first_collection_ref};

my $tol = 0.3;
my @dummy_cluster = ();
my $temp = 2;
my ($first_cluster_ref, $first_substrate_ref) = cluster_iteration(\@first_collection, $tol);
my @first_cluster = @{$first_cluster_ref};
my @first_substrate = @{$first_substrate_ref};

#------------------------------------------------------------
# Using rough approximation to identify the normal direction 

my @atom1 = (0, 0, 0, 0, 0);
my @atom2 = (0, 0, 0, 0, 0);

my $direction = 0;
my $orientation = '0';

for (my $i = 1; $i <= 3; $i++){
	@first_cluster = reverse sort { $a->[$i] <=> $b->[$i] } @first_cluster;
	for (my $j = 0; $j <= 4; $j++){
	    $atom1[$j] = $first_cluster[0][$j];
		$atom2[$j] = $first_cluster[$#first_cluster][$j];
    }
	my ($top_check, $bottom_check) = identify_direction(\@atom1, \@atom2, \@first_substrate, $i);
	if ($top_check eq 'true'){
		$direction = $i;
		$orientation = 'top';
	} 

	if ($bottom_check eq 'true'){
		$direction = $i;
		$orientation = 'bottom';
	}
}

my @true_first_cluster  = @{separate_cluster(\@first_collection, $direction, $orientation, $length)};
@first_substrate = @{identify_substrate(\@first_collection, \@true_first_cluster)};
my $fh;
open($fh, '>', "$ARGV[0]/substrate.in") or die $!;
for (my $i = 0; $i <= $#lattice_vectors; $i++){
	print $fh "$lattice_vectors[$i][0] ","$lattice_vectors[$i][1] ",
	"$lattice_vectors[$i][2] ","$lattice_vectors[$i][3]\n";
}
write_cluster(\@first_substrate, $fh);
close($fh);
open($fh, '>', "$ARGV[0]/cluster.in") or die $!;
write_cluster(\@true_first_cluster, $fh);
close($fh);

