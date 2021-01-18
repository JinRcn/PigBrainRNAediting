#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# To identify uneven RNA editing sites
# Creator: Jinrong Huang <huangjinrong@genomics.cn>
# Date: Nov 1,2019
#
# Three categories are as follows.
# Tissue/Region enriched: At least Delta higher in a particular tissue/region than any other tissue/region
# Group enriched: At least Delta higher in a group of tissues/regions (e.g. 2-4) than any other tissue/region
# Tissue/Region enhanced: At least Delta higher in a particular tissue than the average levels of all tissues/regions
#
# Note: Delta is the difference of editing level, e.g. 0.2 (20%)

my ($in,$ot,$stat,$delta,$rangeGroup,$cutoff);

GetOptions(
	"input=s" => \$in,	# The file including all sites
	"output=s" => \$ot,	# The file including uneven sites
	"stat=s" => \$stat,	# The file summarizing the number of uneven sites
	"delta=s" => \$delta,	# The difference of editing level (range from 0 to 1)
	"rangeGroup=s" => \$rangeGroup,	# The number range for group enriched analysis (e.g. 2-4)
	"cutoff=s" => \$cutoff,	# The minimum level required for an uneven site (e.g. 0.25)
);

$delta ||="0.2";
$rangeGroup ||="2-4";
$cutoff ||="0.25";

my ($min,$max);
if ($rangeGroup =~/(\d+)-(\d+)/){
	$min=$1;
	$max=$2;
}
else {die "Format of rangeGroup should be number-number!\n";}

my %sample;
#my ($TissueEnrichedNumber,$GroupEnrichedNumber,$TissueEnhancedNumber,$other)=("0","0","0","0");
my ($RegionEnrichedNumber,$GroupEnrichedNumber,$RegionEnhancedNumber,$other)=("0","0","0","0");
if($in=~/\.gz$/){open IN,"gunzip -cd <$in|" or die $!;}else{open IN,"<$in" or die $!;}
if($ot=~/\.gz$/){open OT,"|gzip >$ot" or die $!;}else{open OT,">$ot" or die $!;}
while (<IN>){
	# Chromosome_ID   Coordinate      Ref_base        AMY     BG      CB      CC      CTX     HPF     HY      MB      OLF     PM      SC      TH
	# 1       7153    T       NA      NA      0.5     NA      NA      NA      NA      NA      NA      NA      NA      NA
	# 1       2229393 A       NA      0.636363636363636       NA      NA      NA      NA      NA      NA      0.6     NA      NA      NA
	chomp;
	my @t=split /\t/;
	if ($.==1){
		foreach my $i(3..@t-1){
			$sample{$i}=$t[$i];
		}
		print OT "$_\tSortedValue\tMean\tClassification\tRegions\n";
		next;
	}
	my %hash;
	my $line;
	foreach my $i(3..@t-1){
		if ($t[$i] eq "NA"){next;}
#		else {$hash{$i}=$t[$i];}
		else {$hash{$i}=sprintf "%.3f",$t[$i];}
	}
	my (@k,@v,$sum);
	my $n=keys %hash;
	foreach my $k(sort {$hash{$a} <=> $hash{$b}} keys %hash){
		push @k,$k;
		push @v,$hash{$k};
		$sum+=$hash{$k};
		$line.="$hash{$k};";
	}
	$line=~s/;$//;
	my $mean=sprintf "%.3f",$sum/$n;

#	my ($TissueEnriched,$GroupEnriched,$TissueEnhanced);
	my ($RegionEnriched,$GroupEnriched,$RegionEnhanced);

	# TissueEnriched/RegionEnriched
	if ($n eq "1"){
#		$TissueEnriched = $sample{$k[-1]} if $v[-1] >= $cutoff;
		$RegionEnriched = $sample{$k[-1]} if $v[-1] >= $cutoff;
	}
	else {
		my $diff=sprintf "%.3f",$v[-1]-$v[-2];
		if ($diff >= $delta and $v[-1] >= $cutoff){
#			$TissueEnriched = $sample{$k[-1]};
			$RegionEnriched = $sample{$k[-1]};
		}
		else {
			# GroupEnriched
			my @r;
			foreach my $i(0..@v-2){
				my $diff2=sprintf "%.3f",$v[$i+1]-$v[$i];
				push @r,$i+1 if ($diff2 >= $delta and $v[$i+1] >= $cutoff);
			}
			my $r = @v-$r[-1] if defined $r[-1];
			if (defined $r and $r >= $min and $r <= $max){
				foreach my $i($r[-1]..@v-1){
					$GroupEnriched.="$sample{$k[$i]};";
				}
			}
			if (defined $GroupEnriched){$GroupEnriched=~s/;$//;}
			else {
				# TissueEnhanced/RegionEnhanced
				my @h;
				foreach my $i(0..@v-1){
					my $diff2=sprintf "%.3f",$v[$i]-$mean;
					push @h,$i if ($diff2 >= $delta and $v[$i] >= $cutoff);
				}
#				$TissueEnhanced="$sample{$k[$h[0]]}" if @h == "1";
				$RegionEnhanced="$sample{$k[$h[0]]}" if @h == "1";
			}
		}
	}
	my ($type,$detail);
#	if (defined $TissueEnriched){
	if (defined $RegionEnriched){
#		$type="TissueEnriched";
		$type="RegionEnriched";
#		$detail=$TissueEnriched;
		$detail=$RegionEnriched;
#		$TissueEnrichedNumber++;
		$RegionEnrichedNumber++;
	}
	elsif (defined $GroupEnriched){
		$type="GroupEnriched";
		$detail=$GroupEnriched;
		$GroupEnrichedNumber++;
	}
#	elsif (defined $TissueEnhanced){
	elsif (defined $RegionEnhanced){
#		$type="TissueEnhanced";
		$type="RegionEnhanced";
#		$detail=$TissueEnhanced;
		$detail=$RegionEnhanced;
#		$TissueEnhancedNumber++;
		$RegionEnhancedNumber++;
	}
	else {
		$other++;
		next;
	}
	print OT "$_\t$line\t$mean\t$type\t$detail\n";
}
close IN;
close OT;


#print "Type\tNumber\nTissueEnrichedNumber\t$TissueEnrichedNumber\nGroupEnrichedNumber\t$GroupEnrichedNumber\nTissueEnhancedNumber\t$TissueEnhancedNumber\nOhter\t$other\n";
print "Type\tNumber\nRegionEnrichedNumber\t$RegionEnrichedNumber\nGroupEnrichedNumber\t$GroupEnrichedNumber\nRegionEnhancedNumber\t$RegionEnhancedNumber\nOther\t$other\n";

my %Classification2Regions;
if($ot=~/\.gz$/){open IN,"gunzip -cd <$ot|" or die $!;}else{open IN,"<$ot" or die $!;}
while (<IN>){
# Chromosome_ID   Coordinate      Ref_base        AMY     BG      CB      CC      CTX     HPF     HY      MB      OLF     PM      SC      TH      SortedValue     Mean    Classification  Regions
# 1       17148   T       NA      0       NA      0       0.5     0.2     NA      NA      0.0909090909090909      NA      0.416666666666667       NA      0;0;0.0909090909090909;0.2;0.416666666666667;0.5        0.201262626262626       GroupEnriched   SC;CTX
	chomp;
	next if $.==1;
	my @t=split /\t/;
	if ($t[-1]=~/;/){
		my @s=split /;/,$t[-1];
		print "$_\tPlease check rangeGroup!\n" if (@s < $min or @s > $max);
		foreach my $s(@s){
			$Classification2Regions{$t[-2]}{$s}++;
		}
	}
	else {
		$Classification2Regions{$t[-2]}{$t[-1]}++;
	}
}
close IN;

#my @Classification=("TissueEnriched","GroupEnriched","TissueEnhanced");
my @Classification=("RegionEnriched","GroupEnriched","RegionEnhanced");
open OT,">$stat" or die $!;
print OT "Classification\tRegions\tNumber\n";
foreach my $c(@Classification){
	foreach my $k(sort keys %{$Classification2Regions{$c}}){
		print OT "$c\t$k\t$Classification2Regions{$c}{$k}\n";
	}
}
close OT;
exit;
