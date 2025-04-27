use strict;
use warnings;

my $followFile="clinical.txt";
my $expFile="normalizeExp.txt";
#my $sampleFile="sample.txt";
my $sampleFile="all";
my $geneFile="inner.txt";
#my $geneFile="all";

my %geneHash=();
if($geneFile ne 'all')
{
	open(RF,"$geneFile") or die $!;
	while(my $line=<RF>)
	{
		chomp($line);
		$line=~s/\s+//g;
		$geneHash{$line}=1;
	}
	close(RF);
}

my %sampleHash=();
if($sampleFile ne 'all')
{
	open(RF,"$sampleFile") or die $!;
	while(my $line=<RF>)
	{
		chomp($line);
		$line=~s/\s+//g;
		$sampleHash{$line}=1;
	}
	close(RF);
}

my %hash=();
open(RF,"$followFile") or die $!;
while(my $line=<RF>)
{
	next if($line=~/^\n/);
	chomp($line);
	my @arr=split(/\t/,$line);
	my $sampleName=shift(@arr);
  $sampleName=uc($sampleName);
	if($.==1)
	{
		$hash{'id'}=join("\t",@arr);
	}else
	{
	  $hash{$sampleName}=join("\t",@arr);
	}
}
close(RF);

my @sampleName=();
my %expHash=();
my @geneListArr=();
open(RF,"$expFile") or die $!;
while(my $line=<RF>)
{
	chomp($line);
	my @arr=split(/\t/,$line);
	if($.==1)
	{
		@sampleName=@arr;
	}
	elsif($.>2)
	{
	  my @zeroArr=split(/\|\|/,$arr[0]);
		#next if($arr[0]!~/lincRNA/);
		my $flag=0;
		if(($geneFile eq 'all') || (exists $geneHash{$zeroArr[0]}))
		{
			$flag=1;
		}
		if($flag==1)
		{
			push(@geneListArr,$zeroArr[0]);
			for(my $i=1;$i<=$#arr;$i=$i+1)       #需要修改，从第2+1列开始，每隔一列为我们需要的列
			{
				my @subArr=split(/\-/,$sampleName[$i]);
				if($subArr[3]=~/^0/)
				{
					my $subName="$subArr[0]-$subArr[1]-$subArr[2]";
					if(exists $hash{$subName})
					{
						${$expHash{$subName}}{$zeroArr[0]}=$arr[$i];
					}
					#else
					#{
					#	print "$sampleName[$i]\n";
					#}
				}
			}
		}
	}
}
close(RF);

open(WF,">tumor_clinical.txt") or die $!;
print WF "id\t" . $hash{'id'} . "\t" . join("\t",@geneListArr) . "\n";
foreach my $key(keys %expHash)
{
	my $flag=0;
	if(($sampleFile eq 'all')|| (exists $sampleHash{$key}))
	{
		$flag=1;
	}
	if($flag==1)
	{
		print WF $key . "\t" . $hash{$key};
		foreach my $gene(@geneListArr)
		{
			print WF "\t" . ${$expHash{$key}}{$gene};
		}
		print WF "\n";
	}
}
close(WF);

