use strict;
use warnings;
use Math::Trig;


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

	my @cluster  = ();
  my $fname = $_[0];
  open(my $fh, '<:encoding(UTF-8)', $fname)
   or die "Couldn't open file!";
  
  while (my $row = <$fh>){
     chomp $row;
     my @atom = split(' ',$row);
     push(@cluster, [@atom]);
  }
  return \@cluster;
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

sub rotate_cluster_along_x {
  my @cluster = @{$_[0]};
  my $angle   = $_[1];

  my @mean = @{calculate_mean(\@cluster)};
  @cluster = @{shift_to_origin(\@cluster, \@mean)};
  my @rotatedcluster = ();
  my @rotatedatom = (0,0,0,0);
  
  for (my $i = 0; $i <= $#cluster; $i++){
    $rotatedatom[0] = $cluster[$i][0];
    $rotatedatom[1] = $cluster[$i][1];
    $rotatedatom[2] = $cluster[$i][2]*cos($angle) - $cluster[$i][3]*sin($angle);
    $rotatedatom[3] = $cluster[$i][2]*sin($angle) + $cluster[$i][3]*cos($angle); 
    push(@rotatedcluster, [@rotatedatom]);
  }

  @rotatedcluster = @{shift_back_to_original_position(\@rotatedcluster, \@mean)};
  return \@rotatedcluster;
} 

sub rotate_cluster_along_y {
  my @cluster = @{$_[0]};
  my $angle   = $_[1];

  my @mean = @{calculate_mean(\@cluster)};
  @cluster = @{shift_to_origin(\@cluster, \@mean)};
  my @rotatedcluster = ();
  my @rotatedatom = (0,0,0,0);
  
  for (my $i = 0; $i <= $#cluster; $i++){
    $rotatedatom[0] = $cluster[$i][0];
    $rotatedatom[1] = $cluster[$i][1]*cos($angle) - $cluster[$i][3]*sin($angle);
    $rotatedatom[2] = $cluster[$i][2];
    $rotatedatom[3] = $cluster[$i][1]*sin($angle) + $cluster[$i][3]*cos($angle); 
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
  push(@types, $cluster[0][0]);
  for(my $i = 1; $i <= $#cluster; $i++){
    chomp($cluster[$i][0]);
    my $v = 'true';
    for (my $j = ($i - 1); $j >= 0; $j--){
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

sub shift_to_origin {
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
  my @shiftedcluster = @cluster;

  for (my $i = 0; $i <= $#cluster; $i++){
    $shiftedcluster[$i][1] = $cluster[$i][1] + $mean[0];
    $shiftedcluster[$i][2] = $cluster[$i][2] + $mean[1];
    $shiftedcluster[$i][3] = $cluster[$i][3] + $mean[2];
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

sub merge_cuts {
  my @basefirstcut  = @{$_[0]};
  my @basesecondcut = @{$_[1]};
  my @meanfirst  = @{$_[2]};
  my @meansecond = @{$_[3]};
  my @crossover = ();

  @basefirstcut  = @{shift_to_origin(\@basefirstcut, \@meanfirst)};
  @basesecondcut = @{shift_to_origin(\@basesecondcut, \@meansecond)};

  @crossover = @{join_arrays(\@basefirstcut, \@basesecondcut)};
  return \@crossover;
}

# Populating the arrays with file data
#-----------------------------------------------------------------------------------------------

my @firstcluster  = @{read_file($ARGV[0])};
my @secondcluster = @{read_file($ARGV[1])};

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

# Making and Shifting Cuts (Taking out essential atoms) 
#-----------------------------------------------------------------------------------------------

my @finalfirstcut  = @{make_upper_cut(\@firstcluster)};
my @finalsecondcut = @{make_lower_cut(\@secondcluster)};
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
  my @xrotatedcluster  = @{rotate_cluster_along_x(\@firstcluster, $firstclusterangles[0])};
  my @xyrotatedcluster = @{rotate_cluster_along_y(\@xrotatedcluster, $firstclusterangles[1])};

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
  my @xrotatedcluster  = @{rotate_cluster_along_x(\@secondcluster, $secondclusterangles[0])};
  my @xyrotatedcluster = @{rotate_cluster_along_y(\@xrotatedcluster, $secondclusterangles[1])};

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

#-----------------------------------------------------------------------------------------------
# Printing out final crossover result

print "Final Crossover\n";
print_cluster(\@crossover);
my $crossoverlength = $#crossover + 1;
print "Length of Crossover: ", "$crossoverlength\n";
#-----------------------------------------------------------------------------------------------