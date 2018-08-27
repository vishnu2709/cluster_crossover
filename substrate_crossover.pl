use strict;
use warnings;
use Math::Trig;

print "----------------------------------\n";
print "    CROSSOVER ON SUBSTRATE        \n";
print "----------------------------------\n";

# Defining subroutines that are required for the isolating substrate and cluster
# --------------------------------------------------------------------------------

sub read_file {
	my $fname = $_[0];
	my @collection = ();

	open (my $fh, '<:encoding(UTF-8)', $fname)
	 or die "Couldn't open file!";

	while (my $row = <$fh>){
		chomp $row;
		my @atom = split(" ",$row);
		push(@collection, [@atom]);
	}
	return \@collection;
}

sub extract_lattice_vectors {
	my @collection = @{$_[0]};
	my @lattice_vectors = ();
	for (my $i = 0;$i <= 2; $i++){
		push(@lattice_vectors, $collection[$i]);
	}
	for (my $i = 0; $i <= 2; $i++){
		shift @collection;
	}
	return (\@lattice_vectors, \@collection);
}

sub get_magnitude_of_lattice_vectors {
	my @lattice_vectors = @{$_[0]};
	my @magnitude = (0, 0, 0);
	for (my $i = 0; $i <= 2; $i++){
		for (my $j = 1; $j <= 3; $j++){
			$magnitude[$i] = $magnitude[$i] + $lattice_vectors[$i][$j]**2;
		}
		$magnitude[$i] = sqrt($magnitude[$i]);
	}
	return \@magnitude;
}

sub convert_frac_to_atom {
	my @collection = @{$_[0]};
	my @lattice_vectors = @{$_[1]};
	my @magnitude = @{get_magnitude_of_lattice_vectors(\@lattice_vectors)};

	for (my $i = 0; $i <= $#collection; $i++){
		$collection[$i][0] = 'atom';
		for (my $j = 0; $j <= $#magnitude; $j++){
			$collection[$i][$j + 1] = $collection[$i][$j + 1]*$magnitude[$j] 
		}
	}
	return \@collection;
}

sub convert_atom_to_frac {
	my @collection = @{$_[0]};
	my @lattice_vectors = @{$_[1]};
	my @magnitude = @{get_magnitude_of_lattice_vectors(\@lattice_vectors)};

	for (my $i = 0; $i <= $#collection; $i++){
		$collection[$i][0] = 'atom_frac';
		for (my $j = 0; $j <= $#magnitude; $j++){
			$collection[$i][$j + 1] = $collection[$i][$j + 1]/$magnitude[$j] 
		}
	}
	return \@collection;
}

sub relative_coordinates {
	my @atom = @{$_[0]};
	my @collection = @{$_[1]};
	my $index = $_[2];
	my @new_atom = (0, 0, 0, 0, 0);
	my @relative_collection = ();
	splice @collection, $index, 1;
	for (my $i = 0; $i <= $#collection ; $i++){
		$new_atom[0] = $collection[$i][0];
		$new_atom[1] = $collection[$i][1] - $atom[1];
		$new_atom[2] = $collection[$i][2] - $atom[2];
		$new_atom[3] = $collection[$i][3] - $atom[3];
		$new_atom[4] = $collection[$i][4];
		push(@relative_collection, [@new_atom]);
	}
	return \@relative_collection;
}

sub nearest_neighbours {
	my @relative_collection = @{$_[0]};
	my @sorted_neighbours = ();
	my @nearest_neighbour_array = ();
	my @new_atom = (0, 0, 0, 0, 0, 0);
	for (my $i = 0; $i <= $#relative_collection; $i++){
		$new_atom[0] = $relative_collection[$i][0];
		$new_atom[1] = $relative_collection[$i][1];
		$new_atom[2] = $relative_collection[$i][2];
		$new_atom[3] = $relative_collection[$i][3];
		$new_atom[4] = $relative_collection[$i][4];
		for (my $j = 1; $j <= 3; $j++){
			$new_atom[5] = $new_atom[5] + $relative_collection[$i][$j]**2;
		}
		$new_atom[5] = sqrt($new_atom[5]);
		push(@sorted_neighbours, [@new_atom]);
	}
	@sorted_neighbours = sort { $a->[5] <=> $b->[5] } @sorted_neighbours;

	for(my $j = 0; $j <= 5; $j++){
		my @arr = shift @sorted_neighbours;
		push(@nearest_neighbour_array, @arr);
	}
	return \@nearest_neighbour_array;
}

sub lattice_check {
	my @relative_collection = @{$_[0]};
	my $direction = $_[1];
	my $check = 'false';
	for (my $j = 0; $j <= $#relative_collection; $j++){
		if (abs($relative_collection[$j][$direction]) < 0.1){
			$check = 'true';
			next;
		}
	}
	return $check;
}

sub identify_cluster {
	my @collection = @{$_[0]};
	my @cluster = ();
	my $check_a = 0;
	my $check_b = 0;
	my $check_c = 0;
	my @atom = [0, 0, 0, 0, 0];

	for (my $i = 0; $i <= $#collection; $i++){
		$atom[0] = $collection[$i][0];
		$atom[1] = $collection[$i][1];
		$atom[2] = $collection[$i][2];
		$atom[3] = $collection[$i][3];
		$atom[4] = $collection[$i][4];
		my @relative_collection = @{relative_coordinates(\@atom, \@collection, $i)};
		my @nearest_neighbour_array  = @{nearest_neighbours(\@relative_collection)};
		$check_a = lattice_check(\@nearest_neighbour_array, 1);
		$check_b = lattice_check(\@nearest_neighbour_array, 2);
		$check_c = lattice_check(\@nearest_neighbour_array, 3);
		if ($check_a eq 'false' or $check_b eq 'false' or $check_c eq 'false'){
			push(@cluster, $collection[$i]);
		}
	}
	return \@cluster;
}

sub identify_substrate {
	my @collection = @{$_[0]};
	my @cluster = @{$_[1]};

	for (my $j = 0; $j <= $#cluster;$j++){
		for(my $i = 0; $i <= $#collection;$i++){
			if ($cluster[$j][1] == $collection[$i][1] 
				and $cluster[$j][2] == $collection[$i][2] 
				and $cluster[$j][3] == $collection[$i][3]){
				splice @collection, $i, 1;
			}
		}
	}

	return \@collection;
}

# Defining subroutines needed for crossover
# --------------------------------------------------------------------------------

sub print_cluster {
	my @cluster = @{$_[0]};
	for(my $i = 0; $i <= $#cluster; $i++){
		print "$cluster[$i][0] ", "$cluster[$i][1] ", 
		"$cluster[$i][2] ", "$cluster[$i][3] ","$cluster[$i][4]\n";
	}
}

sub calculate_mean {
   my @cluster = @{$_[0]};
   my @mean = (0, 0, 0);

   for(my $i = 0; $i <= $#cluster; $i++){
     for (my $j = 0; $j <= 2; $j++){
        $mean[$j] = $mean[$j] + ($cluster[$i][$j + 1]/($#cluster + 1));
     }
   }
   return \@mean;
}


sub generate_random_angles {
  my @randangles = (0.0, 0.0);
  @randangles[0] = rand(2*pi);
  @randangles[1] = rand(2*pi);
  return \@randangles;
}

sub rotate_cluster_along_a {
  my @cluster = @{$_[0]};
  my @mean = @{$_[1]};
  my $angle   = $_[2];

  @cluster = @{shift_to_new_position(\@cluster, \@mean)};
  my @rotatedcluster = ();
  my @rotatedatom = (0,0,0,0,0);
  
  for (my $i = 0; $i <= $#cluster; $i++){
    $rotatedatom[0] = $cluster[$i][0];
    $rotatedatom[1] = $cluster[$i][1];
    $rotatedatom[2] = $cluster[$i][2]*cos($angle) - $cluster[$i][3]*sin($angle);
    $rotatedatom[3] = $cluster[$i][2]*sin($angle) + $cluster[$i][3]*cos($angle);
    $rotatedatom[4] = $cluster[$i][4]; 
    push(@rotatedcluster, [@rotatedatom]);
  }

  @rotatedcluster = @{shift_back_to_original_position(\@rotatedcluster, \@mean)};
  return \@rotatedcluster;
} 

sub rotate_cluster_along_b {
  my @cluster = @{$_[0]};
  my @mean = @{$_[1]};
  my $angle   = $_[2];

  @mean = @{calculate_mean(\@cluster)};
  @cluster = @{shift_to_new_position(\@cluster, \@mean)};
  my @rotatedcluster = ();
  my @rotatedatom = (0,0,0,0,0);
  
  for (my $i = 0; $i <= $#cluster; $i++){
    $rotatedatom[0] = $cluster[$i][0];
    $rotatedatom[1] = $cluster[$i][1]*cos($angle) - $cluster[$i][3]*sin($angle);
    $rotatedatom[2] = $cluster[$i][2];
    $rotatedatom[3] = $cluster[$i][1]*sin($angle) + $cluster[$i][3]*cos($angle);
    $rotatedatom[4] = $cluster[$i][4];
    push(@rotatedcluster, [@rotatedatom]);
  }
  @rotatedcluster = @{shift_back_to_original_position(\@rotatedcluster, \@mean)};
  return \@rotatedcluster;
}

sub rotate_cluster_along_c {
  my @cluster = @{$_[0]};
  my $angle   = $_[1];

  my @mean = @{calculate_mean(\@cluster)};
  @cluster = @{shift_to_new_position(\@cluster, \@mean)};
  my @rotatedcluster = ();
  my @rotatedatom = (0,0,0,0,0);
  
  for (my $i = 0; $i <= $#cluster; $i++){
    $rotatedatom[0] = $cluster[$i][0];
    $rotatedatom[1] = $cluster[$i][1]*cos($angle) - $cluster[$i][2]*sin($angle);
    $rotatedatom[2] = $cluster[$i][1]*sin($angle) + $cluster[$i][2]*cos($angle);
    $rotatedatom[3] = $cluster[$i][3];
    $rotatedatom[4] = $cluster[$i][4];
    push(@rotatedcluster, [@rotatedatom]);
  }
  @rotatedcluster = @{shift_back_to_original_position(\@rotatedcluster, \@mean)};
  return \@rotatedcluster;
}

sub check_stoichiometry {
  my @cluster = @{$_[0]};
  my @properstoich = @{$_[1]};
  my @types = @{$_[2]};
  my @clusterstoich = @{atoms_of_each_type(\@types, \@cluster)};

  my $output = 'true';
  for (my $i = 0; $i <= $#clusterstoich; $i++){
    if ($clusterstoich[$i] != $properstoich[$i]){
      $output = 'false';
    }
  }
  return $output;
}

sub number_of_types {
  my @types= ();
  my @cluster = @{$_[0]};
  push(@types, $cluster[0][4]);
  for(my $i = 1; $i <= $#cluster; $i++){
    my $v = 'true';
    for (my $j = ($i - 1); $j >= 0; $j--){
      if ($cluster[$i][4] eq $cluster[$j][4]){
        $v = 'false';
      }
    }
    if ($v eq 'true'){
      push(@types, $cluster[$i][4]);
    }
  }
  return \@types;
}

sub atoms_of_each_type {
  my @types = @{$_[0]};
  my @cluster = @{$_[1]};
  my @numberofeachtype = ();
  my $number = 0;

  for (my $i = 0; $i <= $#types; $i++){
    $number = 0;
    for (my $j = 0; $j <= $#cluster; $j++){
      if ($cluster[$j][4] eq $types[$i]){
        $number = $number + 1;
      }
    }
    push(@numberofeachtype, $number);
  }
  return \@numberofeachtype;
}

sub make_upper_cut {
  my @cluster = @{$_[0]};
  my @mean = @{calculate_mean(\@cluster)};
  my @cut  = ();

  for(my $i = 0; $i <= $#cluster; $i++){
    if ($cluster[$i][3] > $mean[2]){
      push(@cut, $cluster[$i]);
    }
  }
  return \@cut;
}

sub make_lower_cut {
  my @cluster = @{$_[0]};
  my @mean = @{calculate_mean(\@cluster)};
  my @cut  = ();

  for(my $i = 0; $i <= $#cluster; $i++){
    if ($cluster[$i][3] < $mean[2]){
      push(@cut, $cluster[$i]);
    }
  }
  return \@cut;
}

sub shift_to_new_position {
  my @cluster = @{$_[0]};
  my @mean = @{$_[1]};
  my @shiftedcluster = @cluster;

  for (my $i = 0; $i <= $#cluster; $i++){
    $shiftedcluster[$i][1] = $cluster[$i][1] - $mean[0];
    $shiftedcluster[$i][2] = $cluster[$i][2] - $mean[1];
    $shiftedcluster[$i][3] = $cluster[$i][3] - $mean[2];
  }
  return \@shiftedcluster;
}

sub shift_back_to_original_position {
  my @cluster = @{$_[0]};
  my @mean = @{$_[1]};
  my @shiftedatom = (0, 0, 0, 0, 0);
  my @shiftedcluster = ();

  for (my $i = 0; $i <= $#cluster; $i++){
  	$shiftedatom[0] = $cluster[$i][0];
    $shiftedatom[1] = $cluster[$i][1] + $mean[0];
    $shiftedatom[2] = $cluster[$i][2] + $mean[1];
    $shiftedatom[3] = $cluster[$i][3] + $mean[2];
    $shiftedatom[4] = $cluster[$i][4];
    push(@shiftedcluster, [@shiftedatom]);
  }
  return \@shiftedcluster;
}

sub join_arrays {
  my @firstarray  = @{$_[0]};
  my @secondarray = @{$_[1]};
  my @crossover = ();
  
  for (my $i = 0; $i <= $#firstarray; $i++){
    push(@crossover, $firstarray[$i]);
  }

  for (my $i = 0; $i <= $#secondarray; $i++){
    push(@crossover, $secondarray[$i]);
  }
  return \@crossover;
}

sub identify_direction {
	my @atom = @{$_[0]};
	my @collection = @{$_[1]};
	my $direction = $_[2];
	my @new_atom = (0, 0, 0, 0, 0);
	my @relative_collection = ();
	for (my $i = 0; $i <= $#collection ; $i++){
		$new_atom[0] = $collection[$i][0];
		$new_atom[1] = $collection[$i][1] - $atom[1];
		$new_atom[2] = $collection[$i][2] - $atom[2];
		$new_atom[3] = $collection[$i][3] - $atom[3];
		$new_atom[4] = $collection[$i][4];
		push(@relative_collection, [@new_atom]);
	}
	my $top_check = 'true';
	my $bottom_check = 'true';
	for (my $i = 0; $i <= $#relative_collection; $i++){
		if ($relative_collection[$i][$direction] > 0){
			$top_check = 'false';
		}
		else {
			$bottom_check = 'false';
		}
	}
	return ($top_check, $bottom_check);
}

sub lowest_value {
	my @cluster = @{$_[0]};
	my $direction = $_[1];
	my $lowest_value = $cluster[0][$direction];
	for (my $i = 1; $i <= $#cluster; $i++){
		if ($cluster[$i][$direction] < $lowest_value){
			$lowest_value = $cluster[$i][$direction];
		}
	}
	return $lowest_value;
}

sub highest_value {
	my @cluster = @{$_[0]};
	my $direction = $_[1];
	my $highest_value = $cluster[0][$direction];
	for (my $i = 1; $i <= $#cluster; $i++){
		if ($cluster[$i][$direction] > $highest_value){
			$highest_value = $cluster[$i][$direction];
		}
	}
	return $highest_value;
}

sub merge_cuts {
  my @basefirstcut  = @{$_[0]};
  my @basesecondcut = @{$_[1]};
  my @meanfirst  = @{$_[2]};
  my @meansecond = @{$_[3]};
  my @crossover  = ();

  my @shiftvector = (0, 0, 0);
  for (my $j = 0; $j <= 2; $j++){
  	$shiftvector[$j] = $meanfirst[$j] - $meansecond[$j];
  }
  @basefirstcut  = @{shift_to_new_position(\@basefirstcut, \@shiftvector)};

  @crossover = @{join_arrays(\@basefirstcut, \@basesecondcut)};
  return \@crossover;
}

# --------------------------------------------------------------------------------
# Reading all of the data into an array

my @first_collection = @{read_file($ARGV[0])};

# Isolating the lattice vectors and finding the magnitude
my ($lattice_vectors_ref, $first_collection_ref)  = extract_lattice_vectors(\@first_collection);
my @lattice_vectors = @{$lattice_vectors_ref};
my @magnitudes = @{get_magnitude_of_lattice_vectors(\@lattice_vectors)};

# Converting fractional coordinates into absolute numbers for simplicity
@first_collection = @{$first_collection_ref};

# Separating cluster from substrate
my @first_cluster   = @{identify_cluster(\@first_collection)};
my @substrate = @{identify_substrate(\@first_collection, \@first_cluster)};

my @second_collection = @{read_file($ARGV[1])};

# Isolating the lattice vectors and finding the magnitude
my ($second_lattice_vectors_ref, $second_collection_ref)  = extract_lattice_vectors(\@second_collection);
@lattice_vectors = @{$second_lattice_vectors_ref};
@magnitudes = @{get_magnitude_of_lattice_vectors(\@lattice_vectors)};

# Converting fractional coordinates into absolute numbers for simplicity
@second_collection = @{$second_collection_ref};

# Separating cluster from substrate
my @second_cluster   = @{identify_cluster(\@second_collection)};

#-----------------------------------------------------------------------------------

print "Original First Cluster\n";
print_cluster(\@first_cluster);
my $firstclusterlength = $#first_cluster + 1; 
print "Length of First Cluster: ", "$firstclusterlength\n\n";

print "Original Second Cluster\n";
print_cluster(\@second_cluster);
my $secondclusterlength = $#second_cluster + 1;
print "Length of Second Cluster: ", "$secondclusterlength\n\n";

#-----------------------------------------------------------------------------------

my @atom = (0, 0, 0, 0, 0);
for (my $j = 0; $j <= $#atom; $j++){
	$atom[$j] = $first_cluster[0][$j];
}

@atom = (0, 0, 0, 0, 0);
my $direction = 0;
my $orientation = '0';

for (my $i = 1; $i <= 3; $i++){
	my ($top_check, $bottom_check) = identify_direction(\@atom, \@substrate, $i);
	if ($top_check eq 'true'){
		$direction = $i;
		$orientation = 'top';
	} 

	if ($bottom_check eq 'true'){
		$direction = $i;
		$orientation = 'bottom';
	}
}
print $direction,",",$orientation,"\n";

my $cluster_lowest_value  = 0;
my $cluster_highest_value = 0;

if ($orientation eq 'top'){
	$cluster_lowest_value = lowest_value(\@first_cluster, $direction);
}
elsif ($orientation eq 'bottom'){
	$cluster_highest_value = highest_value(\@first_cluster, $direction);
}

# Performing Crossover
#------------------------------------------------------------------------------

my @meanfirst = @{calculate_mean(\@first_cluster)};
my @meansecond = @{calculate_mean(\@second_cluster)};

print "Centre Coordinates of First Cluster:\n";
print "$meanfirst[0], ", "$meanfirst[1], ", "$meanfirst[2]\n";
print "Centre Coordinates of Second Cluster:\n";
print "$meansecond[0], ", "$meansecond[1], ", "$meansecond[2]\n\n";

# Finding the types of atoms and no of each
#------------------------------------------------------------------------------

my @types = @{number_of_types(\@first_cluster)};
my @numberofeachtype = @{atoms_of_each_type(\@types, \@first_cluster)};

# Making and Shifting Cuts (Taking out essential atoms) 
#-----------------------------------------------------------------------------------------------

my @finalfirstcut  = @{make_upper_cut(\@first_cluster)};
my @finalsecondcut = @{make_lower_cut(\@second_cluster)};
my @crossover = @{merge_cuts(\@finalfirstcut, \@finalsecondcut, \@meanfirst, \@meansecond)};
my $check     = check_stoichiometry(\@crossover, \@numberofeachtype, \@types);
my $iteration = 0;

my @firstclusterangles   = (0, 0);
my @secondclusterangles  = (0, 0); 
my $finalfirstcutlength  = 0;
my $finalsecondcutlength = 0;

while ($check eq 'false'){
    $iteration = $iteration + 1;
    print "------------------------------\n";
    print "Iteration: ","$iteration\n";
    print "------------------------------\n";
    @firstclusterangles  = @{generate_random_angles()};
    @secondclusterangles = @{generate_random_angles()};

    # First Cluster 
    my @xrotatedcluster  = @{rotate_cluster_along_a(\@first_cluster, 
    	\@meanfirst, $firstclusterangles[0])};
    my @xyrotatedcluster = @{rotate_cluster_along_b(\@xrotatedcluster, 
    	\@meanfirst, $firstclusterangles[1])};

    print "Rotated First Cluster\n";
    print_cluster(\@xyrotatedcluster);
    $firstclusterlength = $#xyrotatedcluster + 1;
    print "\nLength of Rotated First Cluster: ", "$firstclusterlength\n\n";
    @finalfirstcut = @{make_upper_cut(\@xyrotatedcluster)};

    print "First Cut\n";
    print_cluster(\@finalfirstcut);
    $finalfirstcutlength = $#finalfirstcut + 1;
    print "Length of First Cut: ", "$finalfirstcutlength\n\n";

    # Second Cluster 
    @xrotatedcluster  = @{rotate_cluster_along_a(\@second_cluster, 
    	\@meansecond, $secondclusterangles[0])};
    @xyrotatedcluster = @{rotate_cluster_along_b(\@xrotatedcluster, 
    	\@meansecond, $secondclusterangles[1])};

    print "Rotated Second Cluster\n";
    print_cluster(\@xyrotatedcluster);
    $secondclusterlength = $#xyrotatedcluster + 1;
    print "\nLength of Rotated Second Cluster: ", "$secondclusterlength\n\n";
    @finalsecondcut = @{make_lower_cut(\@xyrotatedcluster)};

    print "Second Cut\n";
    print_cluster(\@finalsecondcut);
    $finalsecondcutlength = $#finalsecondcut + 1;
    print "Length of Second Cut: ", "$finalsecondcutlength\n\n";

    @crossover = @{merge_cuts(\@finalfirstcut, \@finalsecondcut, \@meanfirst, \@meansecond)};
    $check     = check_stoichiometry(\@crossover, \@numberofeachtype, \@types);
}

print "Pre-Shift Crossover\n";
print_cluster(\@crossover);
my $crossoverlength = $#crossover + 1;
print "Length of Crossover: ", "$crossoverlength\n";

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

print "Shift Vector\n";
print "$shift_vector[0] ","$shift_vector[1] ","$shift_vector[2]\n";
@crossover = @{shift_to_new_position(\@crossover, \@shift_vector)};

#-----------------------------------------------------------------------------------------------
# Printing out final crossover result

print "Final Crossover\n";
for (my $i = 0; $i <= $#lattice_vectors; $i++){
	print "$lattice_vectors[$i][0] ","$lattice_vectors[$i][1] ",
	"$lattice_vectors[$i][2] ","$lattice_vectors[$i][3]\n";
}
my @new_collection = @{join_arrays(\@substrate, \@crossover)};
# @new_collection = @{convert_atom_to_frac(\@new_collection, \@lattice_vectors)};
print_cluster(\@new_collection);
$crossoverlength = $#crossover + 1;
print "Length of Crossover: ", "$crossoverlength\n";
#-----------------------------------------------------------------------------------------------