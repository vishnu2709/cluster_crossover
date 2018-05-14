use strict;
use warnings;


print " ------------------------------------ \n";
print "          CROSSOVER PROGRAM           \n";
print " ------------------------------------ \n\n";

# Defining subroutines (functions) that are required
#-----------------------------------------------------------------------------------------------

sub print_cluster {
	my @cluster = @{$_[0]};
	for(my $i = 0; $i <= $#cluster; $i++){
		print "atom ", "$cluster[$i][1] ", "$cluster[$i][2] ", "$cluster[$i][3] ","$cluster[$i][0]\n";
	}
}

sub read_file {

	my @firstcluster  = ();
    my @secondcluster = ();

    my $fname = $_[0];
    open(my $fh, '<:encoding(UTF-8)', $fname)
     or die "Couldn't open file!";
    
    while (my $row = <$fh>){
       chomp $row;
       my @atom = split(' ',$row);
       push(@firstcluster, [@atom]);
    }
    
    $fname = $_[1];
    open($fh, '<:encoding(UTF-8)' , $fname)
     or die "Couldn't open file!";
    
    while (my $row = <$fh>){
       chomp $row;
       my @atom = split(' ',$row);
       push(@secondcluster, [@atom]);  
    }
    return (\@firstcluster, \@secondcluster);
}

sub calculate_mean {
   my @firstcluster = @{$_[0]};
   my @secondcluster = @{$_[1]};
   my @meanfirst  = (0, 0, 0);
   my @meansecond = (0, 0, 0);

   for(my $i = 0; $i <= $#firstcluster; $i++){
     for (my $j = 0; $j <= 2; $j++){
        $meanfirst[$j] = $meanfirst[$j] + ($firstcluster[$i][$j + 1]/($#firstcluster + 1));
     }
   }

   for(my $i = 0; $i <= $#secondcluster; $i++){
     for(my $j = 0; $j <= 2; $j++){
        $meansecond[$j] = $meansecond[$j] + ($secondcluster[$i][$j + 1]/($#secondcluster + 1));
     }
   }
   return (\@meanfirst, \@meansecond)
}

sub extract_atoms_of_type {
	my @cluster = @{$_[0]};
	my @types = @{$_[1]};
	my $numbers = @{$_[2]};
	my $count = 0;
	my @extractedcluster = ();
    
    for(my $j = 0; $j <= $#types; $j++){
       for (my $i = 0; $i <= $#cluster; $i++){
       	  if ($type eq $cluster[$i][0]){
       	  	push(@extractedcluster, $cluster[$i]);
       	  	$count = $count + 1
       	  }
       	  if ($count >= $numbers[j]){
       	  	break;
       	  }
       }
   }
}

sub make_cut {
  
  my @firstcluster  = @{$_[0]};
  my @secondcluster = @{$_[1]};
  my $meanfirst  = $_[2];
  my $meansecond = $_[3];    
  my @firstcut   = ();
  my @secondcut  = ();
  
  for(my $i = 0; $i <= $#firstcluster; $i++){
     if ($firstcluster[$i][3] > $meanfirst){
       push(@firstcut, $firstcluster[$i]);
     }
  }

  for(my $i = 0; $i <= $#secondcluster; $i++){
     if ($secondcluster[$i][3] < $meansecond){
       push(@secondcut, $secondcluster[$i]);
     }
   }
   return (\@firstcut, \@secondcut)
}

sub number_of_types {
	my @types= ();
	my @cluster = @{$_[0]};
	push(@types, $cluster[0][0]);
	for(my $i = 1; $i <= $#cluster; $i++){
		chomp($cluster[$i][0]);
		my $v = 'true';
		for (my $j = ($i - 1); $j >= 0; $j--){
            chomp($cluster[$j][0]);
			if ($cluster[$i][0] eq $cluster[$j][0]){
				$v = 'false';
			}
		}
		if ($v eq 'true'){
			push(@types, $cluster[$i][0]);
		}
	}
	return \@types;
}

sub atoms_of_type {
	my @types = @{$_[0]};
	my @cluster = @_{$_[1]};
	my @numberofeachtype = ();

	for (my $i = 0; $i <= $#types; $i++){
		my $number = 0;
		for (my $j = 0; $j ,= $#cluster, $j++){
			if ($cluster[$j][0] eq $types[$i]){
				$number = $number + 1;
			}
		}
		push(@numberofeachtype, $number);
	}
	return \@numberofeachtype;
}

# Defining the arrays
#-----------------------------------------------------------------------------------------------

my @crossover = ();

# Populating the arrays with file data
#-----------------------------------------------------------------------------------------------

my ($firstclusterref, $secondclusterref) = read_file($ARGV[0], $ARGV[1]);
my @firstcluster  = @{$firstclusterref};
my @secondcluster = @{$secondclusterref};

# Now let us calculate the central co-ordinates
#-----------------------------------------------------------------------------------------------

my ($firstmeanref, $secondmeanref) = calculate_mean(\@firstcluster, \@secondcluster);
my @meanfirst = @{$firstmeanref};
my @meansecond = @{$secondmeanref};

print "Centre Coordinates of First Cluster:\n";
print "$meanfirst[0], ", "$meanfirst[1], ", "$meanfirst[2]\n";
print "Centre Coordinates of Second Cluster:\n";
print "$meansecond[0], ", "$meansecond[1], ", "$meansecond[2]\n\n";

# Making cuts (We need to add all the conditions here)
#-----------------------------------------------------------------------------------------------

# 1. Need to first know how many types of atoms are there.
# 2. Then the no of each type has to be calculated.

my @types = @{number_of_types(\@firstcluster)};

my ($firstcutref, $secondcutref) = make_cut(\@firstcluster, \@secondcluster, $meanfirst[2], $meansecond[2]);
my @firstbasecut  = @{$firstcutref};
my @secondbasecut = @{$secondcutref};

#-----------------------------------------------------------------------------------------------
# Printing out cuts

print "Original First Cluster\n";
print_cluster(\@firstcluster);
my $firstclusterlength = $#firstcluster + 1; 
print "Length of First Cluster: ", "$firstclusterlength\n\n";

print "Original Second Cluster\n";
print_cluster(\@secondcluster);
my $secondclusterlength = $#secondcluster + 1;
print "Length of Second Cluster: ", "$secondclusterlength\n\n";

print "First Cut\n";
print_cluster(\@firstcut);
my $firstcutlength = $#firstcut + 1;
print "Length of First Cut: ", "$firstcutlength\n\n";

print "Second Cut\n";
print_cluster(\@secondcut);
my $secondcutlength = $#secondcut + 1;
print "Length of Second Cut: ", "$secondcutlength\n\n";

# Shifting the cuts
#-----------------------------------------------------------------------------------------------

for(my $i = 0; $i <= $#firstcut; $i++){
   for(my $j = 1; $j <= 3; $j++){
      $firstcut[$i][$j] = $firstcut[$i][$j] - $meanfirst[$j - 1];
   }
   push(@crossover, $firstcut[$i]);
}

for(my $i = 0; $i <= $#secondcut; $i++){
   for(my $j = 1; $j <= 3; $j++){
      $secondcut[$i][$j] = $secondcut[$i][$j] - $meansecond[$j - 1];
   }
   push(@crossover, $secondcut[$i]);
}

my $crossoverlength = $#crossover + 1;
print "Length of Crossover, ", "$crossoverlength\n";
print "Final Crossover Result\n";
for(my $i = 0; $i <= $#crossover; $i++){
    print "atom ","$crossover[$i][1] ","$crossover[$i][2] ","$crossover[$i][3] ","$crossover[$i][0]\n";
}
#-----------------------------------------------------------------------------------------------
