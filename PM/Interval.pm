package PM::Interval;

use strict;

my $pkg = __PACKAGE__;
my $minFrac = 0.5; # the minimal region/binSize fraction for a region considered to be valid
my $sliceInterval = 0;
my %includedSeqs;

my $ucsc = '/home/ec2-user/work/data/genomes/ucsc';

my %build2File = (
	'hg19'   =>  "$ucsc/hg19.chrom.sizes",
	'mm10'   =>  "$ucsc/mm10.chrom.sizes"
);

=head2 new

 Title   : new
 Usage   : my $handler=PM::Interval->new(%params);
 Function: create a new interval handler
 Returns : an object of PM::Interval
 Args    : see below

=over 4

=item size
 the size for each interval bin

=item min-frac
 intervals with length/bin-size smaller than this value are skipped in analysis

=item seqs
 an array reference, only consider intervals on these sequences. Default is all sequences.

=item slice-interval
 should intervals provided in a file is sliced into small ones of given size by the option 'size'.
 
=item file
 a file containing all the intervals on which each bin-sized interval will be generated.

=item build
 a genome build name such as 'hg19' or a file containing sequence lengths; this has lower priority than the option 'file'.

=back

=cut

sub new
{
	my ($class,%params) = @_;
	my $self = {};
	bless $self, $pkg;
	if($params{'min-frac'})
	{
		$minFrac = $params{'min-frac'};
	}
	$self->{'_bin_size'} = $params{'size'} if($params{'size'});
	%includedSeqs = map { $_ => 1 } @{$params{'seqs'}} if($params{'seqs'});

	if(my $file = $params{'file'})
	{
		my $fh = _open_file($file) or return undef;
		$self->{'_fh'} = $fh;
		if($params{'slice-interval'}) # break intervals to bin-size
		{
			$sliceInterval = 1;
		}
	}elsif(my $sp = $params{'build'})
	{
		# get chromosome lengths
		my $chrLen = _get_chr_info($sp);
		$self->{'_chr_length'} = $chrLen;
		unless($params{'size'})
		{
			warn "option 'size' is needed for $pkg to create across-genome intervals\n";
			return undef;
		}
	}else
	{
		warn "Can't create an object of $pkg with "._hash_to_str(%params)."\n";
		return undef;
	}
	return $self;
}

=head2 next

 Title   : next
 Usage   : my $i = $handler->next;
 Function: Get next interval
 Returns : an array reference
 Args    : none

=cut

sub _next_child_interval
{
	my $self = shift;
	my $parent = $self->{'_parent_interval'} or 
	die "No parental interval for finding child intervals:$!";
	my ($pChr, $pStart, $pEnd) = @$parent;
	my $binSize = $self->bin_size();
	my $chr;
	my $start;
	my $end;
	if(my $lastInterval = $self->{'_last_interval'}) # parent was visited
	{
		my ($lastChr, $lastStart, $lastEnd) = @$lastInterval;
		if($lastEnd + $binSize > $pEnd) # exceeding range
		{
			my $expectedSize = $pEnd - $lastEnd;
			return () if($expectedSize/$binSize < $minFrac); # too small
			return ($pChr, $lastEnd+1, $pEnd);
		}else
		{
			return ($pChr, $lastEnd+1, $lastEnd + $binSize);
		}

	}else # a new parent
	{
		if($pStart + $binSize -1 > $pEnd)
		{
			my $expectedSize = $pEnd - $pStart + 1;
			return () if($expectedSize/$binSize < $minFrac); # too small
			return ($pChr, $pStart, $pEnd); # the whole region
		}else
		{
			return ($pChr, $pStart, $pStart+$binSize-1);
		}
	}
	
}

sub next
{
	my $self = shift;
	my @interval;
	if(my $fh = $self->{'_fh'})
	{
		# check if a child interval exists
		if($self->{'_parent_interval'})
		{
			@interval = $self->_next_child_interval();
			if(@interval)
			{
				$self->{'_last_interval'} = \@interval;
				return \@interval;
			}
		}
		# otherwise read a new parental interval
		my @fields;
		while(<$fh>) # read a data line
		{
			next if /^#/ or /^\s*$/; # comment or empty line
			chomp;
			@fields = split "\t";
			# check if this sequence is in the included seqs if
			# specified
			next if(%includedSeqs and !$includedSeqs{$fields[0]});
			next unless($fields[1]=~/^\d+$/); 
			die "Line '".$_."' in input interval file doesn't conform to the format ".
			"\nseqName\tstart\tend\n" unless($fields[2]=~/^\d+$/);
			if($sliceInterval) # slice the interval to bin size
			{
				$self->{'_last_interval'} = undef;
				$self->{'_parent_interval'} = [@fields[0..2]];
				@interval = $self->_next_child_interval();
				# this parent interval may be too small to hold any
				# child bin, then
				next unless(@interval);
			}else# otherwise return as a whole
			{
				@interval = @fields[0..2];
			}
			last; # find and get out
		}
		return undef unless(@interval); # no more valid intervals in this file
	}else # extract intervals from a whole genome
	{
		# get the end of last interval first
		my ($chr,$start, $end) = _last_interval($self);
		return undef unless(defined $chr); # no more intervals
		my $chrLen = $self->_get_chr_len($chr);
		my $binSize = $self->bin_size();
		if($end + $binSize > $chrLen) # reach chromosome end
		{
			my $size = $chrLen - $end;
			# test whether the remained region is too short
			if($size/$binSize < $minFrac)
			{
				# get a region from next chromosomes
				while($chr = $self->_next_chr)
				{
					$chrLen = $self->_get_chr_len($chr);
					next if($chrLen/$binSize < $minFrac); # very short chromosome
					last;
				}
				return undef unless(defined $chr); # no more chromosomes
				push @interval, $chr, 1, $chrLen<$binSize? $chrLen : $binSize;

			}else # OK
			{
				push @interval, $chr, $end+1, $chrLen;
			}
		}else
		{
			push @interval, $chr, $end+1, $end+$binSize;
		}
	}
	$self->{'_last_interval'} = \@interval;
	return \@interval;
}

sub _last_interval
{
	my $self = shift;

	return @{$self->{'_last_interval'}} if(exists $self->{'_last_interval'});
	# otherwise the first region in new chromosome
	my $chr = $self->_next_chr or return ();
	return ($chr, undef, 0);
}

sub bin_size
{
	my $self = shift;
	return $self->{'_bin_size'};
}

# internal methods

sub _get_chr_len
{
	my ($self,$chr) = @_;

	return $self->{'_chr_length'}->{$chr};
}

sub _next_chr
{
	my $self = shift;

	unless(exists $self->{'_remaining_chrs'})
	{
		$self->{'_remaining_chrs'} = [sort keys(%{$self->{'_chr_length'}})];
	}

	while(my $chr = shift @{$self->{'_remaining_chrs'}})
	{
		next if(%includedSeqs and !$includedSeqs{$chr});
		return $chr;
	}
	return undef;
}

sub _open_file
{
	my $file = shift;
	my $fh;

	open($fh, "< $file") or die "Can't open $file:$!";

	return $fh;
}

sub _get_chr_info
{
	my $build = shift;

	my $infoFile;
	if(-f $build)
	{
		$infoFile = $build;
	}else
	{
		$infoFile = $build2File{$build} or
		die "The build '$build' has not been implemented yet:$!";
	}

	my $fh = _open_file($infoFile);
	my %chrLens;
	while(<$fh>)
	{
		chomp;
		my ($chr, $len) = split "\t";
		$chrLens{$chr} = $len;
	}
	return \%chrLens;
}

sub _hash_to_str
{
	my %hash = @_;
	my $str = '{';
	while(my ($k, $v) = each %hash)
	{
		$str .= sprintf("%s => %s, ", $k, $v);
	}
	$str .= '}';
	return $str;
}

sub DESTROY
{
	my $self = shift;
	if($self->{'_fh'})
	{
		close $self->{'_fh'};
	}

	$self->SUPER::DESTROY();
}

1;

__END__

