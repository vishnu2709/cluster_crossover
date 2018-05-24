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

sub find_random_distance_along_a {
	my @atom = @{$_[0]};
	my @collection = @{$_[1]};
	my $distance = 0.0;
	for (my $i = 0; $i <= $#collection; $i++){
		if ($atom[2] == $collection[$i][2] and $atom[3] == $collection[$i][3]){
			$distance = $collection[$i][1] - $atom[1];
			return $distance,
		}
	}
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
@collection = @{convert_frac_to_atom(\@collection, \@lattice_vectors)};

for (my $i = 0; $i <= $#lattice_vectors; $i++){
	print "$lattice_vectors[$i][0] ","$lattice_vectors[$i][1] ","$lattice_vectors[$i][2] ","$lattice_vectors[$i][3]\n";
}

for (my $i = 0; $i <= $#collection; $i++){
	print "atom ","$collection[$i][1] ","$collection[$i][2] ","$collection[$i][3] ","$collection[$i][4]\n";
}
