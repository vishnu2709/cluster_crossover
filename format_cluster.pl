use strict;
use warnings;

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

sub print_cluster {
	my @cluster = @{$_[0]};
	for(my $i = 0; $i <= $#cluster; $i++){
		print "atom ", "$cluster[$i][1] ", "$cluster[$i][2] ", "$cluster[$i][3] ","$cluster[$i][0]\n";
	}
}

print_cluster(read_file($ARGV[0]));