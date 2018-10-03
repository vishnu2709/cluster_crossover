#!/usr/bin/perl
use strict;
use warnings;
use Math::Trig;
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use lib dirname(abs_path $0);

require 'crossover_subroutines.pl';
# --------------------------------------------------------------------------------
# Reading all of the data into an array

my $fh;

my $filename = 'substrate_crossover_output.txt';
open($fh, '>', $filename) or die $!;

print $fh "----------------------------------\n";
print $fh "    CROSSOVER ON SUBSTRATE        \n";
print $fh "----------------------------------\n";

my @first_collection = @{read_file($ARGV[0])};

# Isolating the lattice vectors and finding the magnitude
my ($lattice_vectors_ref, $first_collection_ref)  = extract_lattice_vectors(\@first_collection);
my @lattice_vectors = @{$lattice_vectors_ref};
my @magnitudes = @{get_magnitude_of_lattice_vectors(\@lattice_vectors)};

# Converting fractional coordinates into absolute numbers for simplicity
@first_collection = @{$first_collection_ref};

# Separating cluster from substrate
my $tol = 0.3;
my @dummy_cluster = ();
my $temp = 2;
my ($first_cluster_ref, $first_substrate_ref) = cluster_iteration(\@first_collection, $tol);
my @first_cluster = @{$first_cluster_ref};
my @first_substrate = @{$first_substrate_ref};


my @second_collection = @{read_file($ARGV[1])};

# Isolating the lattice vectors and finding the magnitude
my ($second_lattice_vectors_ref, $second_collection_ref)  = extract_lattice_vectors(\@second_collection);
@lattice_vectors = @{$second_lattice_vectors_ref};
@magnitudes = @{get_magnitude_of_lattice_vectors(\@lattice_vectors)};

# Converting fractional coordinates into absolute numbers for simplicity
@second_collection = @{$second_collection_ref};
$tol = 0.3;
# Separating cluster from substrate
my ($second_cluster_ref, $second_substrate_ref) = cluster_iteration(\@second_collection, $tol);
my @second_cluster = @{$second_cluster_ref};
my @second_substrate = @{$second_substrate_ref};

my $length = $ARGV[2];
#-----------------------------------------------------------------------------------

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

print $fh $direction,",",$orientation,"\n";

my @true_first_cluster  = @{separate_cluster(\@first_collection, $direction, $orientation, $length)};
my @true_second_cluster = @{separate_cluster(\@second_collection, $direction, $orientation, $length)};
@first_substrate = @{identify_substrate(\@first_collection, \@true_first_cluster)};

my $cluster_lowest_value  = 0;
my $cluster_highest_value = 0;

if ($orientation eq 'top'){
	$cluster_lowest_value = lowest_value(\@true_first_cluster, $direction);
}
elsif ($orientation eq 'bottom'){
	$cluster_highest_value = highest_value(\@true_first_cluster, $direction);
}


#-----------------------------------------------------------------------------------
print $fh "Original First Cluster\n";
write_cluster(\@true_first_cluster, $fh);
my $firstclusterlength = $#true_first_cluster + 1; 
print $fh "Length of First Cluster: ", "$firstclusterlength\n\n";

print $fh "Original Second Cluster\n";
write_cluster(\@true_second_cluster, $fh);
my $secondclusterlength = $#true_second_cluster + 1;
print $fh "Length of Second Cluster: ", "$secondclusterlength\n\n";

#-----------------------------------------------------------------------------------



# Performing Crossover
#------------------------------------------------------------------------------

my @meanfirst = @{calculate_mean(\@true_first_cluster)};
my @meansecond = @{calculate_mean(\@true_second_cluster)};

print $fh "Centre Coordinates of First Cluster:\n";
print $fh "$meanfirst[0], ", "$meanfirst[1], ", "$meanfirst[2]\n";
print $fh "Centre Coordinates of Second Cluster:\n";
print $fh "$meansecond[0], ", "$meansecond[1], ", "$meansecond[2]\n\n";

# Finding the types of atoms and no of each
#------------------------------------------------------------------------------

my @types = @{number_of_types(\@true_first_cluster)};
my @numberofeachtype = @{atoms_of_each_type(\@types, \@true_first_cluster)};

# Making and Shifting Cuts (Taking out essential atoms) 
#-----------------------------------------------------------------------------------------------

my @finalfirstcut  = @{make_upper_cut(\@true_first_cluster)};
my @finalsecondcut = @{make_lower_cut(\@true_second_cluster)};
my @crossover;
my $check     = 'false';
my $iteration = 0;

my @firstclusterangles   = (0, 0);
my @secondclusterangles  = (0, 0); 
my $finalfirstcutlength  = 0;
my $finalsecondcutlength = 0;

while ($check eq 'false'){
    $iteration = $iteration + 1;
    print $fh "------------------------------\n";
    print $fh "Iteration: ","$iteration\n";
    print $fh "------------------------------\n";
    @firstclusterangles  = @{generate_random_angles()};
    @secondclusterangles = @{generate_random_angles()};

    # First Cluster
    my @xyrotatedcluster = @{generate_rotated_cluster(\@true_first_cluster, 
    	\@meanfirst, \@firstclusterangles, $direction)};

    print $fh "Rotated First Cluster\n";
    write_cluster(\@xyrotatedcluster, $fh);
    $firstclusterlength = $#xyrotatedcluster + 1;
    print $fh "\nLength of Rotated First Cluster: ", "$firstclusterlength\n\n";
    @finalfirstcut = @{make_upper_cut(\@xyrotatedcluster)};

    print $fh "First Cut\n";
    write_cluster(\@finalfirstcut, $fh);
    $finalfirstcutlength = $#finalfirstcut + 1;
    print $fh "Length of First Cut: ", "$finalfirstcutlength\n\n";

    # Second Cluster 
    @xyrotatedcluster = @{generate_rotated_cluster(\@true_second_cluster, 
    	\@meansecond, \@secondclusterangles, $direction)};

    print $fh "Rotated Second Cluster\n";
    write_cluster(\@xyrotatedcluster, $fh);
    $secondclusterlength = $#xyrotatedcluster + 1;
    print $fh "\nLength of Rotated Second Cluster: ", "$secondclusterlength\n\n";
    @finalsecondcut = @{make_lower_cut(\@xyrotatedcluster)};

    print $fh "Second Cut\n";
    write_cluster(\@finalsecondcut, $fh);
    $finalsecondcutlength = $#finalsecondcut + 1;
    print $fh "Length of Second Cut: ", "$finalsecondcutlength\n\n";

    @crossover = @{merge_cuts(\@finalfirstcut, \@finalsecondcut, \@meanfirst, \@meansecond)};
    $check     = check_stoichiometry(\@crossover, \@numberofeachtype, \@types);
}

print $fh "Pre-Shift Crossover\n";
write_cluster(\@crossover, $fh);
my $crossoverlength = $#crossover + 1;
print $fh "Length of Crossover: ", "$crossoverlength\n";

my @shift_vector = (0, 0, 0);
my $new_cluster_lowest_value = 0;
my $new_cluster_highest_value = 0;

if ($orientation eq 'top'){
	$new_cluster_lowest_value = lowest_value(\@crossover, $direction);
    if ($new_cluster_lowest_value < $cluster_lowest_value){
    	$shift_vector[$direction - 1] = $new_cluster_lowest_value - $cluster_lowest_value;
    }
}
elsif ($orientation eq 'bottom'){
	$new_cluster_highest_value = highest_value(\@crossover, $direction);
	if ($new_cluster_highest_value > $cluster_highest_value){
		$shift_vector[$direction - 1] = $new_cluster_highest_value - $cluster_highest_value;
	}
}

print $fh "Shift Vector\n";
print $fh "$shift_vector[0] ","$shift_vector[1] ","$shift_vector[2]\n";
@crossover = @{shift_to_new_position(\@crossover, \@shift_vector)};

#-----------------------------------------------------------------------------------------------
# Printing out final crossover result

print $fh "Final Crossover\n";
for (my $i = 0; $i <= $#lattice_vectors; $i++){
	print "$lattice_vectors[$i][0] ","$lattice_vectors[$i][1] ",
	"$lattice_vectors[$i][2] ","$lattice_vectors[$i][3]\n";
	print $fh "$lattice_vectors[$i][0] ","$lattice_vectors[$i][1] ",
	"$lattice_vectors[$i][2] ","$lattice_vectors[$i][3]\n";
}
my @new_collection = @{join_arrays(\@first_substrate, \@crossover)};
# @new_collection = @{convert_atom_to_frac(\@new_collection, \@lattice_vectors)};
print_cluster(\@new_collection);
write_cluster(\@new_collection, $fh);
$crossoverlength = $#crossover + 1;
print $fh "Length of Crossover: ", "$crossoverlength\n";

close($fh);
#-----------------------------------------------------------------------------------------------
