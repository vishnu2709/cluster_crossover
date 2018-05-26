use strict;
use warnings;
use Math::Trig;

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

sub check_along_a {
	my @atom = @{$_[0]};
	my @collection = @{$_[1]};
	my $check = 'true';
	for (my $j = 0; $j <= $#collection; $j++){
		if ($atom[2] == $collection[$j][2] and $atom[3] == $collection[$j][3]
			and $atom[1] != $collection[$j][1]){
			$check = 'false';
		}
	}
	return $check;
}

sub check_along_b {
	my @atom = @{$_[0]};
	my @collection = @{$_[1]};
	my $check = 'true';
	for (my $j = 0; $j <= $#collection; $j++){
		if ($atom[1] == $collection[$j][1] and $atom[3] == $collection[$j][3]
			and $atom[2] != $collection[$j][2]){
			$check = 'false';
		}
	}
	return $check;
}

sub check_along_c {
	my @atom = @{$_[0]};
	my @collection = @{$_[1]};
	my $check = 'true';
	for (my $j = 0; $j <= $#collection; $j++){
		if ($atom[1] == $collection[$j][1] and $atom[2] == $collection[$j][2]
			and $atom[3] != $collection[$j][3]){
			$check = 'false';
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
	my @atom = (0, 0, 0, 0, 0);

	for (my $i = 0; $i <= $#collection; $i++){
		$atom[0] = $collection[$i][0];
		$atom[1] = $collection[$i][1];
		$atom[2] = $collection[$i][2];
		$atom[3] = $collection[$i][3];
		$atom[4] = $collection[$i][4];
		$check_a = check_along_a(\@atom, \@collection);
		$check_b = check_along_b(\@atom, \@collection);
		$check_c = check_along_c(\@atom, \@collection);
		if ($check_a eq 'true' and $check_b eq 'true' and $check_c eq 'true'){
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

@collection = @{$collection_ref};

my @magnitude = @{get_magnitude_of_lattice_vectors(\@lattice_vectors)};

my @cluster = @{separate_cluster_from_substrate(\@collection)};

for (my $i = 0; $i <= $#cluster; $i++){
	print "$cluster[$i][0] ","$cluster[$i][1] ","$cluster[$i][2] ","$cluster[$i][3] ","$cluster[$i][4]\n"; 
}
print "$#cluster\n";
print "$#collection\n";