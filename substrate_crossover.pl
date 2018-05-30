use strict;
use warnings;
use Math::Trig;
use POSIX qw(ceil floor);

print "----------------------------------\n";
print "    CROSSOVER ON SUBSTRATE        \n";
print "----------------------------------\n";

# Defining subroutines that are required for the crossover
# ---------------------------------------------------------

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

sub separate_cluster_from_substrate {
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
		$check_a = lattice_check(\@relative_collection, 1);
		$check_b = lattice_check(\@relative_collection, 2);
		$check_c = lattice_check(\@relative_collection, 3);
		if ($check_a eq 'false' or $check_b eq 'false' or $check_c eq 'false'){
			push(@cluster, $collection[$i]);
		}
	}
	return \@cluster;
}

sub convert_frac_to_atom {
	my @collection = @{$_[0]};
	my @lattice_vectors = @{$_[1]};
	my @magnitude = @{get_magnitude_of_lattice_vectors(\@lattice_vectors)};

	for (my $i = 0; $i <= $#collection; $i++){
		for (my $j = 0; $j <= $#magnitude; $j++){
			$collection[$i][$j + 1] = $collection[$i][$j + 1]*$magnitude[$j] 
		}
	}
	return \@collection;
}

# ---------------------------------------------------------

my @collection = @{read_file("geometry.in")};

my ($lattice_vectors_ref, $collection_ref)  = extract_lattice_vectors(\@collection);
my @lattice_vectors = @{$lattice_vectors_ref};
my @magnitudes = @{get_magnitude_of_lattice_vectors(\@lattice_vectors)};

@collection = @{convert_frac_to_atom($collection_ref, $lattice_vectors_ref)};

my @cluster = @{separate_cluster_from_substrate(\@collection)};
print "Third Pass\n\n";
for (my $i = 0; $i <= $#cluster; $i++){
	print "$cluster[$i][0] ","$cluster[$i][1] ","$cluster[$i][2] ","$cluster[$i][3] ","$cluster[$i][4]\n"; 
}
print "$#cluster\n";
print "$#collection\n";