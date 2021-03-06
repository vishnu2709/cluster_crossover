#!/usr/bin/perl
use strict;
use warnings;
use Math::Trig;

# Defining subroutines that are required for the isolating substrate and cluster
# --------------------------------------------------------------------------------

sub read_file {
	my $fname = $_[0];
	my @collection = ();

	open GEOMETRY,  "< $fname" or die $!;

	while (my $row = <GEOMETRY>){
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
	my $direction_1 = $_[1];
	my $direction_2 = $_[2];
	my $tolerance = $_[3];
	my $check = 'false';
	for (my $j = 0; $j <= $#relative_collection; $j++){
		if (abs($relative_collection[$j][$direction_1]) < $tolerance and 
			abs($relative_collection[$j][$direction_2]) < $tolerance){
			$check = 'true';
			next;
		}
	}
	return $check;
}

sub identify_cluster {
	my @collection = @{$_[0]};
	my $tolerance = $_[1];
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
		$check_a = lattice_check(\@relative_collection, 1, 2, $tolerance);
		$check_b = lattice_check(\@relative_collection, 2, 3, $tolerance);
		$check_c = lattice_check(\@relative_collection, 1, 3, $tolerance);
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

sub write_cluster {
	my @cluster = @{$_[0]};
	my $fhandle = $_[1];
	for(my $i = 0; $i <= $#cluster; $i++){
		print $fhandle "$cluster[$i][0] ", "$cluster[$i][1] ",
		"$cluster[$i][2] ", "$cluster[$i][3] ", "$cluster[$i][4]\n";
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
  $randangles[0] = rand(2*pi);
  $randangles[1] = rand(2*pi);
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

sub generate_rotated_cluster {
	my @cluster   = @{$_[0]};
	my @mean      = @{$_[1]};
	my @angles    = @{$_[2]};
	my $direction = $_[3];
	my @xrotatedcluster = ();
	my @xyrotatedcluster = ();

	if ($direction == 1) {
		@xrotatedcluster  = @{rotate_cluster_along_b(\@cluster, 
    	\@mean, $angles[0])};
        @xyrotatedcluster = @{rotate_cluster_along_c(\@xrotatedcluster, 
    	\@mean, $angles[1])};
	}
	elsif ($direction == 2) {
        @xrotatedcluster  = @{rotate_cluster_along_a(\@cluster, 
    	\@mean, $angles[0])};
        @xyrotatedcluster = @{rotate_cluster_along_c(\@xrotatedcluster, 
    	\@mean, $angles[1])};
	}
	else {
		@xrotatedcluster  = @{rotate_cluster_along_a(\@cluster, 
    	\@mean, $angles[0])};
        @xyrotatedcluster = @{rotate_cluster_along_b(\@xrotatedcluster, 
    	\@mean, $angles[1])};
	}
	return \@xyrotatedcluster;
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
  my @shiftedatom = (0, 0, 0, 0, 0);
  my @shiftedcluster = ();

  for (my $i = 0; $i <= $#cluster; $i++){
  	$shiftedatom[0] = $cluster[$i][0];
    $shiftedatom[1] = $cluster[$i][1] - $mean[0];
    $shiftedatom[2] = $cluster[$i][2] - $mean[1];
    $shiftedatom[3] = $cluster[$i][3] - $mean[2];
    $shiftedatom[4] = $cluster[$i][4];
    push(@shiftedcluster, [@shiftedatom]);
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

sub separate_cluster {
	my @collection = @{$_[0]};
	my $direction = $_[1];
	my $orientation = $_[2];
	my $length = $_[3];
	my @cluster = ();
	my @atom = (0, 0, 0, 0, 0);
	if ($orientation eq 'top') {
		@collection = reverse sort { $a->[$direction] <=> $b->[$direction] } @collection;
		for (my $i = 0; $i < $length; $i++){
			$atom[0] = $collection[$i][0];
			$atom[1] = $collection[$i][1];
			$atom[2] = $collection[$i][2];
			$atom[3] = $collection[$i][3];
			$atom[4] = $collection[$i][4];
			push(@cluster, [@atom]);
		}
	}
	elsif ($orientation eq 'bottom') {
		@collection = sort { $a->[$direction] <=> $b->[$direction] } @collection;
		for (my $i = 0; $i < $length; $i++){
			$atom[0] = $collection[$i][0];
			$atom[1] = $collection[$i][1];
			$atom[2] = $collection[$i][2];
			$atom[3] = $collection[$i][3];
			$atom[4] = $collection[$i][4];
			push(@cluster, [@atom]);
		}
	}
	return \@cluster;
}

sub identify_direction {
	my @atom1 = @{$_[0]};
	my @atom2 = @{$_[1]};
	my @collection = @{$_[2]};
	my $direction = $_[3];
	my @new_atom = (0, 0, 0, 0, 0);
	my @relative_collection1 = ();
	my @relative_collection2 = ();
	for (my $i = 0; $i <= $#collection ; $i++){
		$new_atom[0] = $collection[$i][0];
		$new_atom[1] = $collection[$i][1] - $atom1[1];
		$new_atom[2] = $collection[$i][2] - $atom1[2];
		$new_atom[3] = $collection[$i][3] - $atom1[3];
		$new_atom[4] = $collection[$i][4];
		push(@relative_collection1, [@new_atom]);
	}
	for (my $i = 0; $i <= $#collection ; $i++){
		$new_atom[0] = $collection[$i][0];
		$new_atom[1] = $collection[$i][1] - $atom2[1];
		$new_atom[2] = $collection[$i][2] - $atom2[2];
		$new_atom[3] = $collection[$i][3] - $atom2[3];
		$new_atom[4] = $collection[$i][4];
		push(@relative_collection2, [@new_atom]);
	}
	my $top_check = 'true';
	my $bottom_check = 'true';
	for (my $i = 0; $i <= $#relative_collection1; $i++){
		if ($relative_collection1[$i][$direction] > 0){
			$top_check = 'false';
		}
	}
	for (my $i = 0; $i <= $#relative_collection2; $i++){
		if ($relative_collection2[$i][$direction] < 0){
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

sub cluster_iteration {
	my @collection = @{$_[0]};
	my $tolerance = $_[1];
	my $cluster_ref = identify_cluster(\@collection, $tolerance);
	my $substrate_ref = identify_substrate(\@collection, $cluster_ref);
	return ($cluster_ref, $substrate_ref);
}

1;