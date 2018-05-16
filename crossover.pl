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

	  my @a  = ();
    my @b = ();

    my $fname = $_[0];
    open(my $fh, '<:encoding(UTF-8)', $fname)
     or die "Couldn't open file!";
    
    while (my $row = <$fh>){
       chomp $row;
       my @atom = split(' ',$row);
       push(@a, [@atom]);
    }
    
    $fname = $_[1];
    open($fh, '<:encoding(UTF-8)' , $fname)
     or die "Couldn't open file!";
    
    while (my $row = <$fh>){
       chomp $row;
       my @atom = split(' ',$row);
       push(@b, [@atom]);  
    }
    return (\@a, \@b);
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

sub generate_random_numbers {
  my @numberofeachtype = @{$_[0]};
  my @randnumbers = ();
  my @complementarynumbers = ();
  my $tempnumber = 0;
  for (my $i = 0; $i <= $#numberofeachtype; $i++){
    $tempnumber = int(rand($numberofeachtype[$i] + 1));
    push(@randnumbers, $tempnumber);
    push(@complementarynumbers, ($numberofeachtype[$i] - $tempnumber))
  }
  return (\@randnumbers, \@complementarynumbers);
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

sub atoms_of_each_type {
  my @types = @{$_[0]};
  my @cluster = @{$_[1]};
  my @numberofeachtype = ();
  my $number = 0;

  for (my $i = 0; $i <= $#types; $i++){
    $number = 0;
    for (my $j = 0; $j <= $#cluster; $j++){
      if ($cluster[$j][0] eq $types[$i]){
        $number = $number + 1;
      }
    }
    push(@numberofeachtype, $number);
  }
  return \@numberofeachtype;
}

sub extract_atoms_of_type {
	my @cluster = @{$_[0]};
	my @types = @{$_[1]};
	my @numbers = @{$_[2]};
	my $count = 0;
	my @extractedcluster = ();
    
    for(my $j = 0; $j <= $#types; $j++){
       $count = 0;
       for (my $i = 0; $i <= $#cluster; $i++){
          if ($count >= $numbers[$j]){
            last;
          }
       	  if ($types[$j] eq $cluster[$i][0]){
       	  	push(@extractedcluster, $cluster[$i]);
       	  	$count = $count + 1;
       	  }
       }
   }
   return \@extractedcluster;
}

sub shift_upwards {
  my @cluster = @{$_[0]};
  my $lowest = 0;

  for (my $i = 0; $i <= $#cluster; $i++){
    if ($cluster[$i][3] < $lowest){
      $lowest = $cluster[$i][3];
    }
  }
  for (my $i = 0; $i <= $#cluster; $i++){
    if ($lowest != 0){
       $cluster[$i][3] = $cluster[$i][3] - $lowest + 0.1;
    }
  }
  return \@cluster;
}

sub shift_downwards {
  my @cluster = @{$_[0]};
  my $highest = 0;

  for (my $i = 0; $i <= $#cluster; $i++){
    if ($cluster[$i][3] > $highest){
      $highest = $cluster[$i][3];
    }
  }
  for (my $i = 0; $i <= $#cluster; $i++){
    if ($highest != 0){
      $cluster[$i][3] = $cluster[$i][3] - $highest - 0.1;
    }
  }
  return \@cluster;
}

sub shift_to_origin {
  my @cluster = @{$_[0]};
  my @mean = @{calculate_mean(\@cluster)};

  for (my $i = 0; $i <= $#cluster; $i++){
    $cluster[$i][1] = $cluster[$i][1] - $mean[1];
    $cluster[$i][2] = $cluster[$i][2] - $mean[2];
  }
  return \@cluster;
}

sub join_cuts {
  my @firstcut = @{$_[0]};
  my @secondcut = @{$_[1]};
  my @crossover = ();
  
  for (my $i = 0; $i <= $#firstcut; $i++){
    push(@crossover, $firstcut[$i]);
  }

  for (my $i = 0; $i <= $#secondcut; $i++){
    push(@crossover, $secondcut[$i]);
  }
  return \@crossover;
}

# Populating the arrays with file data
#-----------------------------------------------------------------------------------------------

my ($firstclusterref, $secondclusterref) = read_file($ARGV[0], $ARGV[1]);
my @firstcluster  = @{$firstclusterref};
my @secondcluster = @{$secondclusterref};
@firstcluster  = reverse sort { $a->[3] <=> $b->[3] } @firstcluster;
@secondcluster = sort { $a->[3] <=> $b->[3] } @secondcluster;

print "Original First Cluster\n";
print_cluster(\@firstcluster);
my $firstclusterlength = $#firstcluster + 1; 
print "Length of First Cluster: ", "$firstclusterlength\n\n";

print "Original Second Cluster\n";
print_cluster(\@secondcluster);
my $secondclusterlength = $#secondcluster + 1;
print "Length of Second Cluster: ", "$secondclusterlength\n\n";

# Now let us calculate the central co-ordinates
#-----------------------------------------------------------------------------------------------

my @meanfirst = @{calculate_mean(\@firstcluster)};
my @meansecond = @{calculate_mean(\@secondcluster)};

print "Centre Coordinates of First Cluster:\n";
print "$meanfirst[0], ", "$meanfirst[1], ", "$meanfirst[2]\n";
print "Centre Coordinates of Second Cluster:\n";
print "$meansecond[0], ", "$meansecond[1], ", "$meansecond[2]\n\n";

# Finding the types of atoms and no of each
#-----------------------------------------------------------------------------------------------

my @types = @{number_of_types(\@firstcluster)};
my @numberofeachtype = @{atoms_of_each_type(\@types, \@firstcluster)};

# Generating a random number of atoms of each type to extract from the cluster
#-----------------------------------------------------------------------------------------------

my ($randomnumbersref, $complementarynumbersref) = generate_random_numbers(\@numberofeachtype);
my @randnumbers = @{$randomnumbersref};
my @complementarynumbers = @{$complementarynumbersref}; 

# Base First Cut (Taking out essential atoms) 
#-----------------------------------------------------------------------------------------------

my @basefirstcut = @{extract_atoms_of_type(\@firstcluster, \@types, \@randnumbers)};
my @basesecondcut = @{extract_atoms_of_type(\@secondcluster, \@types, \@complementarynumbers)};

# Shifting the cuts and joining them to form the final crossover result
#-----------------------------------------------------------------------------------------------

my @upwardshiftedfirstcut    = @{shift_upwards(\@basefirstcut)};
my @downwardshiftedsecondcut = @{shift_downwards(\@basesecondcut)};

my @finalfirstcut = @{shift_to_origin(\@upwardshiftedfirstcut)};
my @finalsecondcut = @{shift_to_origin(\@downwardshiftedsecondcut)};

my @crossover = @{join_cuts(\@finalfirstcut, \@finalsecondcut)};

#-----------------------------------------------------------------------------------------------
# Printing out cuts

print "First Cut\n";
print_cluster(\@finalfirstcut);
my $finalfirstcutlength = $#finalfirstcut + 1;
print "Length of First Cut: ", "$finalfirstcutlength\n\n";

print "Second Cut\n";
print_cluster(\@finalsecondcut);
my $finalsecondcutlength = $#finalsecondcut + 1;
print "Length of Second Cut: ", "$finalsecondcutlength\n\n";

print "Final Crossover\n";
print_cluster(\@crossover);
my $crossoverlength = $#crossover + 1;
print "Length of Crossover: ", "$crossoverlength\n";
#-----------------------------------------------------------------------------------------------
